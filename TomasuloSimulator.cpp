#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define NUM_CONFIGS 9
#define MAX_CONFIG_LEN 19
#define MAX_COMMAND_LEN 8
#define MAX_CONFIG_FILE_SIZE 1024
#define MEMIN_INSTRUCTIONS_NUM 4096
#define MAX_MEMIN_FILE_SIZE (10*MEMIN_INSTRUCTIONS_NUM)
#define MAX_FILE_LINE_LEN 50
#define NUM_REGISTERS 16
#define ADD_ID 0
#define MUL_ID 1
#define DIV_ID 2
#define ADD_OPCODE 2
#define SUB_OPCODE 3
#define MUL_OPCODE 4
#define DIV_OPCODE 5
#define HALT_OPCODE 6 
#define INSTRUCTION_QUEUE_SIZE 16

typedef struct {
    uint32_t cdb_id;        // Identifier for the CDB (ADD, MUL, DIV)
    uint32_t output;        // output transmitted on the CDB
    uint32_t busy;          // Flag indicating if the station is busy (1) or available (0)
    uint32_t cycle;         // The cycle in which the data is being delivered
    uint32_t pc;            // The program counter of the instruction whose execution resulted in writing to the CDB
} cdb_t;

typedef struct {
    uint32_t cdb_id;        // Identifier for the CDB (e.g., ADD, MUL, DIV)
    uint32_t tag;           // number of the CDB in the group
    uint32_t busy;          // Flag indicating if the station is busy (1) or available (0)
    uint32_t pc;            // The program counter of instruction the unit is executing
    uint32_t delay;         // Cost
} unit_t;

typedef struct {
    uint32_t station_id;    // Identifier for the reservation station
    uint32_t busy;          // Flag indicating if the station is busy (1) or available (0)
    uint32_t ins;           // instruction number
    uint32_t op;            // Operation to be performed (e.g., ADD, MUL, DIV)
    uint32_t vj;            // Value of operand j
    uint32_t vk;            // Value of operand k
    uint32_t qj;            // Reservation station producing Vj
    uint32_t qk;            // Reservation station producing Vk
    uint32_t start_cycle;   // The cycle in which the current instruction has started to execute. Relevant only if busy=1. If value is 0, the instruction has not started yet.
} reservation_station_t;

typedef struct {
    long raw_instruction;       // The instruction as it appears in the memin file
    uint32_t opcode;            // Opcode of the instruction
    uint32_t dst;               // Destination register index
    uint32_t src0;              // Source register 0 index
    uint32_t src1;              // Source register 1 index
    uint32_t pc;                // The program counter of the instruction
    uint32_t tag;               // The tag of the station the instruction was executed in
    uint32_t cycle_issued;      // Cycle when the instruction was issued
    uint32_t cycle_execute_start; // Cycle when execution of the instruction started
    uint32_t cycle_execute_end; // Cycle when execution of the instruction ended
    uint32_t cycle_cdb;         // Cycle when the result was written to the Common Data 
    reservation_station_t *reservation_station;     // Pointer to the reservation station the instruction is waiting in
                                                    // NULL if unassigned
} instruction_t;

typedef struct configurations
{ //why long?
    long add_nr_units;
    long mul_nr_units;
    long div_nr_units;
    long add_nr_reservation;
    long mul_nr_reservation;
    long div_nr_reservation;
    long add_delay;
    long mul_delay;
    long div_delay;
} configurations_t;

typedef struct {
    float* reg; // Array to store 16 single-precision floating-point registers
    uint32_t pc;   // Program Counter register
    reservation_station_t* add_reservation_stations; // Dynamically allocated array of reservation stations
    reservation_station_t* mul_reservation_stations; // Dynamically allocated array of reservation stations
    reservation_station_t* div_reservation_stations; // Dynamically allocated array of reservation stations
    cdb_t* cdb_add; // add Data Buses
    cdb_t* cdb_mul; // mul Data Buses
    cdb_t* cdb_div; // div Data Buses
    uint32_t cycle;
    configurations_t conf;
    unit_t* add_units; //  Dynamically allocated array of units
    unit_t* mul_units; //  Dynamically allocated array of units
    unit_t* div_units; //  Dynamically allocated array of units
    instruction_t* instruction_que; // Array to store up to 16 instructions 
    uint8_t instruction_que_size;   // Indicates the current instruction queue size.
    uint8_t halted;                 // Value of 1 means that the halt command was executed. Value of 0 otherwise.
} processor_t;

typedef struct {
    float v_i;      // The stable value of the reqister
    uint32_t q_i;   // The id of the station that the register is waiting for
    uint8_t busy;   // 0 if value is valid, 1 if going to be "2" next cycle, 2 if waiting for a reservation station
} register_t;

// Configuation names as they appear in the configuration input file
const char CONFIGS[NUM_CONFIGS][MAX_CONFIG_LEN] = {
    "add_nr_units",
    "mul_nr_units",
    "div_nr_units",
    "add_nr_reservation",
    "mul_nr_reservation",
    "div_nr_reservation",
    "add_delay",
    "mul_delay",
    "div_delay"
};

void init_configs(configurations_t* configurations)
{
    configurations->add_nr_units = -1;
    configurations->mul_nr_units = -1;
    configurations->div_nr_units = -1;
    configurations->add_nr_reservation = -1;
    configurations->mul_nr_reservation = -1;
    configurations->div_nr_reservation = -1;
    configurations->add_delay = -1;
    configurations->mul_delay = -1;
    configurations->div_delay = -1;
}

void init_CDB(cdb_t* cdb, uint32_t cdb_id, uint32_t tag)        // output transmitted on the CDB
{
    //
}

void init_instruction(instruction_t* o_instruction)
{
    // TODO init all fields
}

void init_processor(processor_t* processor, configurations_t* configs, uint32_t* memin_data, size_t memin_size) {
    // Initialize processor state variables
    
    // Initialize registers based on configuration data
    // TODO init correctly according to register_t
    // TODO why use malloc?
    processor->reg = (float*)malloc(sizeof(float) * NUM_REGISTERS);
    if (processor->reg == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize registers to zero or other default values
    for (int i = 0; i < 16; i++)
    {
        processor->reg[i].v_i = float(i);
    }
    // Initialize reservation stations based on configuration data
    int num_reservation_stations = configs->add_nr_reservation + configs->mul_nr_reservation + configs->div_nr_reservation;
    processor->reservation_stations = (reservation_station_t*)malloc(sizeof(reservation_station_t) * num_reservation_stations);
    if (processor->reservation_stations == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize reservation stations
    for (int i = 0; i < num_reservation_stations; i++)
    {
        // TODO first station_id should be 1, not 0.
    }
    // Initialize units based on configuration data
    int num_units = configs->add_nr_units + configs->mul_nr_units + configs->div_nr_units;
    processor->units = (unit_t*)malloc(sizeof(unit_t) * num_units);
    if (processor->units == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize CDBs
    for (int i = 0; i < num_reservation_stations; i++)
    {
        //
    }
    // Copy data from memin_data to memory locations in the processor

    // Additional initialization based on parsed data

    // Free memory allocated for configurations after initialization
    // TODO what memory was allocated??
    free(configs);
}

void init_reservation_station(reservation_station_t* reservation_station)
{
    // TODO
}

/*
 * Function: print_configs
 * ------------------------
 * Prints the parsed configuration settings to the console.
 *
 * configurations: Pointer to the configurations structure containing the parsed settings.
 */
void print_configs(const configurations_t configurations)
{
    printf("Configs:\n");
    printf("\tadd_nr_units = %ld\n", configurations.add_nr_units);
    printf("\tmul_nr_units = %ld\n", configurations.mul_nr_units);
    printf("\tdiv_nr_units = %ld\n", configurations.div_nr_units);
    printf("\tadd_nr_reservation = %ld\n", configurations.add_nr_reservation);
    printf("\tmul_nr_reservation = %ld\n", configurations.mul_nr_reservation);
    printf("\tdiv_nr_reservation = %ld\n", configurations.div_nr_reservation);
    printf("\tadd_delay = %ld\n", configurations.add_delay);
    printf("\tmul_delay = %ld\n", configurations.mul_delay);
    printf("\tdiv_delay = %ld\n", configurations.div_delay);
    printf("\n");
}

/*
 * Function: parse_configs
 * -----------------------
 * Parses a configuration file to extract and set configuration settings for the processor.
 *
 * config_file: Pointer to the configuration file to be parsed.
 * configurations: Pointer to the configurations structure where parsed settings will be stored.
 *
 * Returns:
 *  - 0: Successful parsing of the configuration file.
 *  - Non-zero: Indicates an error during memory allocation, file reading, or parsing.
 */
int parse_configs(FILE* config_file, configurations_t* configurations)
{
    char* file_contents_allocated_buffer = (char*)malloc(MAX_CONFIG_FILE_SIZE);
    char* file_contents = file_contents_allocated_buffer;
    if (file_contents == NULL)
    {
        printf("Failed to allocate!\n");
        return 1;
    }
    memset(file_contents, 0, MAX_CONFIG_FILE_SIZE);

    char line[MAX_FILE_LINE_LEN] = { 0 };

    long value = 0;
    size_t line_len = 0;

    size_t file_size = fread(file_contents, 1, MAX_CONFIG_FILE_SIZE, config_file);
    // TODO check for errors?
    if (file_size == MAX_CONFIG_FILE_SIZE)
    {
        printf("File too big!\n");
        return -1;
    }

    size_t content_index = 0;
    while (1)
    {
        content_index = 0;
        while (file_contents[content_index] != '\n')
        {
            ++content_index;
        }
        // include the newline character
        ++content_index;

        memcpy(line, file_contents, content_index);

        for (size_t i = 0; i < NUM_CONFIGS; i++)
        {
            // find which config is the correct one
            if (strstr(line, CONFIGS[i]) != NULL)
            {
                // found the config, now parse the value for this config
                // 3 is the len of " = "
                value = strtol(line + strnlen(CONFIGS[i], MAX_CONFIG_LEN) + 3, NULL, 10);
                // NOTE: not checking conversion errors! if an error occurred here, return 1

                // assign the input value to the configuration struct
                switch (i)
                {
                case 0:
                    configurations->add_nr_units = value;
                    break;
                case 1:
                    configurations->mul_nr_units = value;
                    break;
                case 2:
                    configurations->div_nr_units = value;
                    break;
                case 3:
                    configurations->add_nr_reservation = value;
                    break;
                case 4:
                    configurations->mul_nr_reservation = value;
                    break;
                case 5:
                    configurations->div_nr_reservation = value;
                    break;
                case 6:
                    configurations->add_delay = value;
                    break;
                case 7:
                    configurations->mul_delay = value;
                    break;
                case 8:
                    configurations->div_delay = value;
                    break;
                default:
                    break;
                }
                // found the config and assigned the value, move to next config
                break;
            }
        }

        memset(line, 0, sizeof(line));
        file_contents += content_index;
        if (file_contents_allocated_buffer + file_size <= file_contents)
        {
            break;
        }
    }

    free(file_contents_allocated_buffer);
    return 0;
}

/*
 * Function: parse_memin
 * ---------------------
 * Parses a memory initialization file to extract and process memory initialization instructions.
 *
 * memin_file: Pointer to the memory initialization file to be parsed.
 *
 * Returns:
 *  - 0: Successful parsing of the memory initialization file.
 *  - Non-zero: Indicates an error during memory allocation, file reading, or parsing.
 */
int parse_memin(FILE* memin_file, instruction_t* instructions)
{
    char* file_contents_allocated_buffer = (char*)malloc(MAX_MEMIN_FILE_SIZE);
    char* file_contents = file_contents_allocated_buffer;

    if (file_contents == NULL)
    {
        printf("Failed to allocate!\n");
        return 1;
    }
    memset(file_contents, 0, MAX_MEMIN_FILE_SIZE);

    char line[MAX_FILE_LINE_LEN] = { 0 };

    long raw_instruction = 0;
    uint8_t opcode = 0;
    uint8_t dst = 0;
    uint8_t src0 = 0;
    uint8_t src1 = 0;

    long opcode_mask = 0x0f000000;
    long dst_mask = 0x00f00000;
    long src0_mask = 0x000f0000;
    long src1_mask = 0x0000f000;

    size_t file_size = fread(file_contents, 1, MAX_MEMIN_FILE_SIZE, memin_file);
    // TODO check for errors?
    if (file_size == MAX_MEMIN_FILE_SIZE)
    {
        printf("File too big!\n");
        return -1;
    }

    uint32_t instruction_index = 0;
    size_t content_index = 0;
    while (1)
    {
        content_index = 0;
        while (file_contents[content_index] != '\n')
        {
            ++content_index;
        }
        // include the newline character
        ++content_index;

        memcpy(line, file_contents, content_index);

        raw_instruction = strtol(line, NULL, 16);
        instructions[instruction_index].pc = instruction_index;
        instructions[instruction_index].raw_instruction = raw_instruction;
        instructions[instruction_index].opcode = (raw_instruction & opcode_mask) >> 24;
        instructions[instruction_index].dst = (raw_instruction & dst_mask) >> 20;
        instructions[instruction_index].src0 = (raw_instruction & src0_mask) >> 16;
        instructions[instruction_index].src1 = (raw_instruction & src1_mask) >> 12;
        ++instruction_index;
        printf("opcode: %u; dst: %u; src0: %u; src1: %u\n", instructions[instruction_index].opcode, instructions[instruction_index].dst, instructions[instruction_index].src0, instructions[instruction_index].src1);

        memset(line, 0, sizeof(line));

        file_contents += content_index;
        if (file_contents_allocated_buffer + file_size <= file_contents)
        {
            break;
        }
    }

    free(file_contents_allocated_buffer);
    return 0;
}

int execute_instructions() {
    // Implementation of executing the instructions based on configurations and memory initialization
    // Return 0 on success, non-zero value on failure
    return 0;
}

int write_regout_file(FILE* regout_file) {
    // Implementation of writing the register output file
    // Return 0 on success, non-zero value on failure
    return 0;

}

int write_traceinst_file(FILE* traceinst_file) {
    // Implementation of writing the trace instruction file
    // Return 0 on success, non-zero value on failure
    return 0;
    return 0;
}

int write_tracecdb_file(FILE* tracecdb_file) {
    // Implementation of writing the trace CDB file
    // Return 0 on success, non-zero value on failure
    return 0;
}

void init_reservation_station_type(reservation_station_t* stations, int num_stations) {
    for (int i = 0; i < num_stations; i++) {
        stations[i].station_id = i;
        stations[i].busy = 0;
        stations[i].ins = 0;
        stations[i].op = 0;
        stations[i].vj = 0;
        stations[i].vk = 0;
        stations[i].qj = 0;
        stations[i].qk = 0;
    }
}

void init_reservation_stations(processor_t* processor, configurations_t* configs) {
    // Initialize add reservation stations
    processor->add_reservation_stations = init_reservation_stations_type(configs->add_nr_reservation);

    // Initialize mul reservation stations
    processor->mul_reservation_stations = init_reservation_stations_type(configs->mul_nr_reservation);

    // Initialize div reservation stations
    processor->div_reservation_stations = init_reservation_stations_type(configs->div_nr_reservation);
}

reservation_station_t* init_reservation_stations_type(int num_stations) {
    reservation_station_t* stations = (reservation_station_t*)malloc(sizeof(reservation_station_t) * num_stations);
    if (stations == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    init_reservation_station_type(stations, num_stations);
    return stations;
}

void reg_init(processor_t *processor)
{
    processor->reg = (float *)malloc(sizeof(float) * NUM_REGISTERS);
    if (processor->reg == NULL)
    {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize registers to zero or other default values
    for (int i = 0; i < 16; i++)
void fetch_instructions(instruction_t* instructions, processor_t* processor)
{
    uint8_t instructions_to_copy = 0;

    // Instruction queue is full - nothing to do
    if (processor->instruction_que_size == INSTRUCTION_QUEUE_SIZE)
    {
        return;
    }

    if (processor->instruction_que_size == INSTRUCTION_QUEUE_SIZE - 1)
    {
        // Only one place in the instruction queue - fetch only one instruction
        instructions_to_copy = 1;
    }
    else
    {
        // At least two places in the instruction queue - fetch two instructions
        instructions_to_copy = 2;
    }

    memcpy(processor->instruction_que + processor->instruction_que_size, instructions, instructions_to_copy);
    processor->instruction_que_size += instructions_to_copy;

}

void start_execution_if_possible(processor_t* processor, uint32_t cycle) {
    // Implementation of executing the instructions based on configurations and memory initialization
    // Return 0 on success, non-zero value on failure

    // TODO where should we write the execute and write cdb cycles?
    // the instruction queue gets renewed all the time, the processor only stores the reservation stations, the memin struct is already out of reach and it seems incorrect to store it there.
    
    // This function must be called on the clock cycle AFTER the issue of the instructions
    uint32_t num_reservation_stations = processor->conf.add_nr_reservation + processor->conf.mul_nr_reservation + processor->conf.div_nr_reservation;
    uint32_t num_execution_units = processor->conf.add_nr_units + processor->conf.mul_nr_units + processor->conf.div_nr_units;
    
    for (uint32_t i = 0; i < num_reservation_stations; ++i)
    {
        // make sure that the station is not busy and not waiting for any other instruction and not already started execution
        if (processor->reservation_stations[i].busy != 0 || 
            processor->reservation_stations[i].qj != 0 || 
            processor->reservation_stations[i].qk != 0 ||
            processor->reservation_stations[i].start_cycle != 0)
        {
            continue;
        }
        
        for (uint32_t j = 0; j < num_execution_units; j++)
        {
            if (processor->units[j].busy == 1)
            {
                continue;
            }
            
            uint32_t opcode = processor->reservation_stations[i].op;
            uint32_t unit_id = processor->units[j].cdb_id;
            if (
                !((opcode == ADD_OPCODE || opcode == SUB_OPCODE) && unit_id == ADD_ID) &&
                !(opcode == MUL_OPCODE && unit_id == MUL_ID) &&
                !(opcode == DIV_OPCODE && unit_id == DIV_ID)
                )
            {
                continue;
            }
            
            processor->units[j].busy = 1;
        }
        
        processor->reservation_stations[i].start_cycle = cycle;
    }
}

uint8_t execute_task_if_finished_and_get_result(uint32_t cycle, reservation_station_t reservation_station, configurations_t conf, float* result)
{
    long delay = 0;
    
    switch (reservation_station.op)
    {
    case ADD_OPCODE:
        *result = reservation_station.vj + reservation_station.vk;
        delay = conf.add_delay;
        break;
    case SUB_OPCODE:
        *result = reservation_station.vj - reservation_station.vk;
        delay = conf.add_delay;
        break;
    case MUL_OPCODE:
        *result = reservation_station.vj * reservation_station.vk;
        delay = conf.mul_delay;
        break;
    case DIV_OPCODE:
        if (reservation_station.vk == 0)
        {
            printf("Division by zero!\n");
            exit(-1);
        }
        *result = reservation_station.vj / reservation_station.vk;
        delay = conf.div_delay;
        break;
    case HALT_OPCODE:
        //TODO
        // *halted = 1; // ???
        printf("\n\nTODO!\n\n");
        exit(-2);
    default:
        exit(11);
        break;
    }

    return (cycle == reservation_station.start_cycle + delay) ? 1 : 0;
}

void write_cdb(processor_t *processor, uint32_t cycle)
{
    uint32_t nr_reservation = 
        processor->conf.add_nr_reservation + 
        processor->conf.mul_nr_reservation + 
        processor->conf.div_nr_reservation;
    uint32_t nr_units =
        processor->conf.add_nr_units +
        processor->conf.mul_nr_units +
        processor->conf.div_nr_units;
    
    float result = 0;
    
    for (uint32_t i = 0; i < nr_reservation; ++i)
    {
        processor->reg[i] = float(i);
    }
}

void init_cdb(cdb_t* cdb, int cdb_id) {
    cdb->cdb_id = cdb_id;
    cdb->busy = 0;
    cdb->output = 0;
}

void cdb_init(processor_t* processor) {
    // Initialize CDBs based on configuration data
    cdb_t* add = (cdb_t*)malloc(sizeof(cdb_t));
    cdb_t* mul = (cdb_t*)malloc(sizeof(cdb_t));
    cdb_t* div = (cdb_t*)malloc(sizeof(cdb_t));
    if (add == NULL || mul == NULL || div == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    init_cdb(add, ADD_ID);
    init_cdb(mul, MUL_ID);
    init_cdb(div, DIV_ID);
    processor->cdb_add = add;
    processor->cdb_mul = mul;
    processor->cdb_div = div;
}

void ALU_unit_init(processor_t *processor, configurations_t *configs)
{
    processor->add_units = (unit_t *)malloc(sizeof(unit_t) * configs->add_nr_units);
    processor->mul_units = (unit_t *)malloc(sizeof(unit_t) * configs->mul_nr_units);
    processor->div_units = (unit_t *)malloc(sizeof(unit_t) * configs->div_nr_units);
    if (processor->add_units == NULL || processor->mul_units == NULL || processor->div_units == NULL)
    {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
}

void instruction_queue_init(processor_t *processor)
{
    for (size_t i = 0; i < 16; i++)
    {
        processor->instruction_que[i].opcode = 0;
        processor->instruction_que[i].dst = 0;
        processor->instruction_que[i].src0 = 0;
        processor->instruction_que[i].src1 = 0;
        processor->instruction_que[i].cycle_issued = 0;
        processor->instruction_que[i].cycle_execute_start = 0;
        processor->instruction_que[i].cycle_execute_end = 0;
        processor->instruction_que[i].cycle_cdb = 0;
    }
}

void init_processor(processor_t *processor, configurations_t *configs, uint32_t *memin_data, size_t memin_size)
{
    // Initialize processor state variables
    processor->cycle = 0;
    processor->pc = 0;
    processor->conf = *configs;
    processor->instruction_que = (instruction_t*)malloc(sizeof(instruction_t) * 16);
    if (processor->instruction_que == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize reservation stations based on configuration data
    init_reservation_stations(processor, configs);
    // Initialize registers based on configuration data
    reg_init(processor);
    // Initialize units based on configuration data
    ALU_unit_init(processor, configs);
    // Initialize CDBs based on configuration data
    cdb_init(processor);
    // initialize the instruction queue
    instruction_queue_init(processor);
    // Copy data from memin_data to memory locations in the processor
    for (size_t i = 0; i < memin_size; i++) {
        processor->instruction_que[i].opcode = (memin_data[i] & 0x0f000000) >> 24;
        processor->instruction_que[i].dst = (memin_data[i] & 0x00f00000) >> 20;
        processor->instruction_que[i].src0 = (memin_data[i] & 0x000f0000) >> 16;
        processor->instruction_que[i].src1 = (memin_data[i] & 0x0000f000) >> 12;
    }
    // Additional initialization based on parsed data
    // Copy data from memin_data to memory locations in the processor
    
    
    // Additional initialization based on parsed data
    
    // Free memory allocated for configurations after initialization
    free(configs);
        if (processor->reservation_stations[i].busy == 0 || processor->reservation_stations[i].start_cycle == 0)
        {
            continue;
        }

        if (execute_task_if_finished_and_get_result(cycle, processor->reservation_stations[i], processor->conf, &result) == 0)
        {
            continue;
        }
        
        // Go over awaiting reservation stations and update with result
        for (uint32_t j = 0; j < nr_reservation; ++j)
        {
            if (i == j)
            {
                continue;
            }

            if (processor->reservation_stations[j].qj == processor->reservation_stations[i].station_id)
            {
                processor->reservation_stations[j].vj = result;
            }

            if (processor->reservation_stations[j].qk == processor->reservation_stations[i].station_id)
            {
                processor->reservation_stations[j].vk = result;
            }
        }

        // Go over awaiting registers and update with result
        for (uint32_t k = 0; k < NUM_REGISTERS; k++)
        {
            if (processor->reg[k].q_i == processor->reservation_stations[i].station_id)
            {
                processor->reg[k].v_i = result;
            }
             
        }

        // TODO where to store this?
        instruction_t instruction;
        instruction.cycle_cdb = cycle;

        // mark execution unit as not busy
        for (uint32_t j = 0; j < nr_units; j++)
        {
            if (processor->reservation_stations[i].ins == processor->units[j].pc)
            {
                processor->units[j].busy = 0;
            }
        }

        // empty reservation station
        processor->reservation_stations[i].busy = 0;
        init_reservation_station(&processor->reservation_stations[i]);
    }

}

// TODO probably irrelevant
void check_finished_tasks(processor_t processor, uint32_t cycle)
{
    uint8_t *add_stations_finished = (uint8_t*)malloc(processor.conf.add_nr_reservation * sizeof(uint8_t));
    uint8_t *mul_stations_finished = (uint8_t*)malloc(processor.conf.mul_nr_reservation * sizeof(uint8_t));
    uint8_t *div_stations_finished = (uint8_t*)malloc(processor.conf.div_nr_reservation * sizeof(uint8_t));

    // create a list of size of all stations and for each maintain a list whether it's done or not, so you can check what stations/registers are waiting for it.
    // don't forget to manage memory correctly
    
    // TODO loop once! write less code by checking the current instrution/type of the reservation station and doing switch-case
    // and assigning the delay (and other distignguishing params) there instead of inside the loop

    // Iterate over add stations
    for (uint32_t i = 0; i < processor.conf.add_nr_reservation; ++i)
    {
        if (processor.reservation_stations[i].busy == 0 || processor.reservation_stations[i].start_cycle == 0)
        {
            continue;
        }
        
        if (cycle == processor.reservation_stations[i].start_cycle + processor.conf.add_delay)
        {
            add_stations_finished[i] = 1;
        }

        if (processor.reservation_stations[i].busy)
        {

        }
    }
    
    // Iterate over mul stations
    for (uint32_t i = processor.conf.add_nr_reservation; i < processor.conf.add_nr_reservation + processor.conf.mul_nr_reservation; ++i)
    {
        if (processor.reservation_stations[i].busy == 0 || processor.reservation_stations[i].start_cycle == 0)
        {
            continue;
        }

        if (cycle == processor.reservation_stations[i].start_cycle + processor.conf.mul_delay)
        {
            mul_stations_finished[i] = 1;
        }

        if (processor.reservation_stations[i].busy)
        {

        }
    }
    
    // Iterate over div stations
    for (uint32_t i = processor.conf.add_nr_reservation + processor.conf.mul_nr_reservation; i < processor.conf.add_nr_reservation + processor.conf.mul_nr_reservation + processor.conf.div_nr_reservation; ++i)
    {
        if (processor.reservation_stations[i].busy == 0 || processor.reservation_stations[i].start_cycle == 0)
        {
            continue;
        }

        if (cycle == processor.reservation_stations[i].start_cycle + processor.conf.div_delay)
        {
            div_stations_finished[i] = 1;
        }

        if (processor.reservation_stations[i].busy)
        {

        }
    }

}

void get_reservation_stations_indices_by_opcode(uint32_t opcode, configurations_t configurations, uint32_t *o_start_index, uint32_t *o_end_index)
{
    switch (opcode)
    {
    case ADD_OPCODE:
    case SUB_OPCODE:
        *o_start_index = 0;
        *o_end_index = *o_start_index + configurations.add_nr_reservation;
        break;
    case MUL_OPCODE:
        *o_start_index = configurations.add_nr_reservation;
        *o_end_index = *o_start_index + configurations.mul_nr_reservation;
        break;
    case DIV_OPCODE:
        *o_start_index = configurations.add_nr_reservation + configurations.mul_nr_reservation;
        *o_end_index = *o_start_index + configurations.div_nr_reservation;
        break;
    case HALT_OPCODE:
        // TODO: what should we do here?
    default:
        printf("try_to_issue_instructions: Got invalid instruction");
        exit(-1);
    }
}

void issue_instructions(processor_t* processor, instruction_t* o_first_instruction, instruction_t* o_second_instruction, const uint32_t current_cycle)
    // Get instrtucions from queue to the reservation stations, if possible, and mark the newly occupied reservation stations as busy.
    // Also, write the station id of reservation station inside the instruction
    // After issuing at most 2 instructions, reorganize the instruction queue
{
    if (processor == NULL)
    {
        printf("try_to_issue_instructions: Got NULL!\n");
        exit(-1);
    }
    
    // if there are no instructions to issue
    if (processor->instruction_que_size == 0)
    {
        return;
    }

    // TODO: This can't happen unless these are allocated in this function
    //o_first_instruction = NULL;
    //o_second_instruction = NULL;

    uint32_t start_point = 0;
    uint32_t end_point = processor->conf.add_nr_reservation;
    uint32_t num_instructions_issued = 0;

    get_reservation_stations_indices_by_opcode(processor->instruction_que[0].opcode, processor->conf, &start_point, &end_point);
    for (uint32_t i = start_point; i < end_point; ++i)
    {
        if (processor->reservation_stations[i].busy == 0)
        {
            *o_first_instruction = processor->instruction_que[0];
            processor->reservation_stations[i].busy = 1;
            o_first_instruction->reservation_station = &processor->reservation_stations[i];
            num_instructions_issued = 1;
            o_first_instruction->cycle_issued = current_cycle;
            break;
        }
    }

    get_reservation_stations_indices_by_opcode(processor->instruction_que[1].opcode, processor->conf, &start_point, &end_point);
    for (uint32_t i = start_point; i < end_point; ++i)
    {
        if (processor->reservation_stations[i].busy == 0)
        {
            *o_second_instruction = processor->instruction_que[1];
            processor->reservation_stations[i].busy = 1;
            o_second_instruction->reservation_station = &processor->reservation_stations[i];
            num_instructions_issued = 2;
            o_second_instruction->cycle_issued = current_cycle;
            break;
        }
    }

    // Update instruction queue by moving the next instructions to the beginning of the queue and updating the correct size
    memmove(processor->instruction_que, processor->instruction_que + num_instructions_issued, INSTRUCTION_QUEUE_SIZE - num_instructions_issued);
    processor->instruction_que_size -= num_instructions_issued;
    memset(processor->instruction_que + processor->instruction_que_size, 0, num_instructions_issued);
}

void assign_register_values_if_possible(register_t* regs, instruction_t instruction)
{
    if (regs[instruction.src0].busy != 2)
    {
        instruction.reservation_station->vj = regs[instruction.src0].v_i;
        instruction.reservation_station->qj = 0;
    }
    else
    {
        instruction.reservation_station->qj = regs[instruction.src0].q_i;
        instruction.reservation_station->vj = UINT16_MAX;
    }

    if (regs[instruction.src1].busy != 2)
    {
        instruction.reservation_station->vk = regs[instruction.src1].v_i;
        instruction.reservation_station->qk = 0;
    }
    else
    {
        instruction.reservation_station->qk = regs[instruction.src1].q_i;
        instruction.reservation_station->vk = UINT16_MAX;
    }

    regs[instruction.dst].busy = 1;
    regs[instruction.dst].q_i = instruction.reservation_station->station_id;
}

void assign_instruction_to_reservation_station(processor_t *processor, instruction_t instruction)
{
    // This instruction didn't get assigned to a reservation station, skip it.
    if (instruction.reservation_station == NULL)
    {
        return;
    }
    
    // this function edits all the parameters of the reservation station to which a new instruction got assigned.
    // it checks whether the required registers are ready, and if it's ready - a next function call will start the execution of the
    // instruction, if possible.

    instruction.reservation_station->ins = instruction.pc;
    instruction.reservation_station->op = instruction.opcode;
    assign_register_values_if_possible(processor->reg, instruction);
}

void set_pending_register_to_busy(register_t *regs)
{
    for (uint32_t i = 0; i < NUM_REGISTERS; ++i)
    {
        regs[i].busy = (regs[i].busy == 0) ? 0 : 2;
    }
}

uint8_t check_if_can_exit(reservation_station_t* reservation_stations)
{
    // TODO

    return 0;
}

int run_processor(processor_t* processor, instruction_t* instructions)
{
    uint32_t cycle = 0;

    instruction_t first_instruction, second_instruction;
    init_instruction(&first_instruction);
    init_instruction(&second_instruction);

    // Fetch first two instructions
    // TODO is this redundant? or should we first fetch and then loop
    fetch_instructions(instructions, processor);

    while (1)
    {
        // Internally, all these functions do their job only when everything is ready for it
        
        // Should we first fetch and then issue? can it happen the other way around?

        // The reservation stations list will be constructed like this:
            // first add_nr_res_stations are the adder stations
            // the mul_nr_res_stations are the multiplier stations
            // the the last div_nr_res_stations are the divider stations

        start_execution_if_possible(processor, cycle);

        write_cdb(processor, cycle);

        //get_instructions_to_execute(processor);
        //execute(processor);
        
        //write_tracecdb_file(cycle);
    
        // if got halted, continue execution loop but don't fetch and issue new instructions.
        // loop until there are no more instructions to issue
        if (processor->halted == 1)
        {
            if (check_if_can_exit(processor->reservation_stations) == 1)
            {
                break;
            }
            continue;
        }

        fetch_instructions(instructions, processor);
        issue_instructions(processor, &first_instruction, &second_instruction, cycle);

        assign_instruction_to_reservation_station(processor, first_instruction);
        assign_instruction_to_reservation_station(processor, second_instruction);

        init_instruction(&first_instruction);
        init_instruction(&second_instruction);

        set_pending_register_to_busy(processor->reg);
        ++cycle;
    }

    //write_traceinst_file(/*instructions*/);
    //write_regout_file(/*processor*/);
}

// TODO probably irrelevant
void execute(processor_t* processor, instruction_t instruction)
{
    const long opcode = instruction.opcode;
    const uint32_t dst = instruction.dst;
    const uint32_t src0 = instruction.src0;
    const uint32_t src1 = instruction.src1;

    switch (opcode)
    {
    case ADD_OPCODE:
        processor->reg[dst] = processor->reg[src0] + processor->reg[src1];
        break;
    case SUB_OPCODE:
        processor->reg[dst] = processor->reg[src0] - processor->reg[src1];
        break;
    case MUL_OPCODE:
        processor->reg[dst] = processor->reg[src0] * processor->reg[src1];
        break;
    case DIV_OPCODE:
        processor->reg[dst] = processor->reg[src0] / processor->reg[src1];
        break;
    case HALT_OPCODE:
        processor->halted = 1;
    default:
        exit(10);
        break;
    }
}



/*
 * Function: main
 * --------------
 * The main function of the program that handles command-line arguments and orchestrates the execution flow.
 *
 * argc: Number of command-line arguments passed to the program.
 * argv: Array of strings containing the command-line arguments.
 *
 * Returns:
 *  - 0: Successful execution of the program.
 *  - Non-zero: Indicates an error during program execution.
 */
int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        printf("Wrong param num %d, should be 3.\n", argc);
        return 1;
    }

    FILE* config_file = NULL;
    FILE* memin_file = NULL;
    // FILE *regout_file = NULL;
    // FILE *traceinst_file = NULL;
    // FILE *tracecdb_file = NULL;
    configurations_t configs;

    errno_t err_config_file = 0;
    errno_t err_memin_file = 0;
    errno_t err_regout_file = 0;
    errno_t err_traceinst_file = 0;
    errno_t err_tracecdb_file = 0;

    instruction_t instructions[MEMIN_INSTRUCTIONS_NUM] = { 0 };

    err_config_file = fopen_s(&config_file, argv[1], "r");
    if (config_file == NULL)
    {
        printf("Failed to open file");
        return 2;
    }

    if (parse_configs(config_file, &configs) != 0)
    {
        printf("Failed to parse configs");
        fclose(config_file);
        return 3;
    }
    fclose(config_file);
    print_configs(configs);


    err_memin_file = fopen_s(&memin_file, argv[2], "r");
    if (memin_file == NULL)
    {
        printf("Failed to open file");
        return 2;
    }

    if (parse_memin(memin_file, instructions) != 0)
    {
        printf("Failed to parse memin");
        fclose(memin_file);
        return 3;
    }
    fclose(memin_file);

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
