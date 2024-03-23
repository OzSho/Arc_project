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

// TODO fix exit codes
// TODO check variable type warnings
// TODO add NO_STATION_ID constant

typedef struct {
    uint32_t cdb_id;        // Identifier for the CDB (ADD, MUL, DIV)
    uint32_t station_id;    // The identifier of the reservation station that sent the data
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
{ 
    unsigned long add_nr_units;
    unsigned long mul_nr_units;
    unsigned long div_nr_units;
    unsigned long add_nr_reservation;
    unsigned long mul_nr_reservation;
    unsigned long div_nr_reservation;
    unsigned long add_delay;
    unsigned long mul_delay;
    unsigned long div_delay;
} configurations_t;

typedef struct {
    float v_i;      // The stable value of the reqister
    uint32_t q_i;   // The id of the station that the register is waiting for
    uint8_t busy;   // 0 if value is valid, 1 if going to be "2" next cycle, 2 if waiting for a reservation station
} register_t;

typedef struct {
    register_t* reg; // Array to store 16 single-precision floating-point registers
    reservation_station_t* reservation_stations; // Dynamically allocated array of reservation stations
    cdb_t cdb_add; // add Data Buses
    cdb_t cdb_mul; // mul Data Buses
    cdb_t cdb_div; // div Data Buses
    uint32_t cycle;
    configurations_t conf;
    unit_t* units; //  Dynamically allocated array of units
    instruction_t* instruction_que; // Array to store up to 16 instructions 
    uint8_t instruction_que_size;   // Indicates the current instruction queue size.
    uint8_t halted;                 // Value of 1 means that the halt command was executed. Value of 0 otherwise.
} processor_t;


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

/**
 * @brief Initializes the configurations structure with default values.
 *
 * @param configurations Pointer to the configurations structure to be initialized.
 */
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

/**
 * @brief Initializes the Common Data Bus (CDB) structure with default values.
 *
 * @param cdb Pointer to the CDB structure to be initialized.
 * @param cdb_id The ID of the CDB.
 */
void init_CDB(cdb_t* cdb, uint32_t cdb_id)
{
    cdb->cdb_id = cdb_id;
    cdb->busy = 0;
    cdb->station_id = 0;
    cdb->output = 0;
    cdb->cycle = 0;
    cdb->pc = 0;
}

/**
 * @brief Initializes an instruction with default values.
 *
 * @param o_instruction Pointer to the instruction to be initialized.
 */
void init_instruction(instruction_t* o_instruction)
{
    o_instruction->raw_instruction = 0;
    o_instruction->opcode = 0;
    o_instruction->dst = 0;
    o_instruction->src0 = 0;
    o_instruction->src1 = 0;
    o_instruction->pc = 0;
    o_instruction->tag = 0;
    o_instruction->cycle_issued = 0;
    o_instruction->cycle_execute_start = 0;
    o_instruction->cycle_execute_end = 0;
    o_instruction->cycle_cdb = 0;
    o_instruction->reservation_station = NULL;
}

/**
 * @brief Initializes a reservation station with default values.
 *
 * @param reservation_station Pointer to the reservation station to be initialized.
 * @param station_id The ID of the reservation station.
 */
void init_reservation_station(reservation_station_t* reservation_station, uint32_t station_id)
{
    reservation_station->station_id = station_id;
    reservation_station->busy = 0;
    reservation_station->ins = 0;
    reservation_station->op = 0;
    reservation_station->vj = 0;
    reservation_station->vk = 0;
    reservation_station->qj = 0;
    reservation_station->qk = 0;
    reservation_station->start_cycle = 0;
}

/**
 * @brief Initializes a unit with default values.
 *
 * @param unit Pointer to the unit to be initialized.
 * @param cdb_id The ID of the Common Data Bus (CDB).
 */
void init_unit(unit_t* unit, uint32_t cdb_id)
{
    unit->cdb_id = cdb_id;
    unit->tag = 0;
    unit->busy = 0;
    unit->pc = 0;
    unit->delay = 0;
}

/**
 * @brief Initializes the processor structure with default values.
 *
 * @param processor Pointer to the processor structure to be initialized.
 * @param configs The configurations for the processor.
 */
void init_processor(processor_t* processor, configurations_t configs/*, uint32_t* memin_data, size_t memin_size*/) {
    // Initialize processor state variables
    
    // TODO no need to use malloc
    // Initialize registers based on configuration data
    processor->reg = (register_t*)malloc(sizeof(register_t) * NUM_REGISTERS);
    if (processor->reg == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize registers to zero or other default values
    for (int i = 0; i < 16; i++)
    {
        processor->reg[i].q_i = 0;
        processor->reg[i].v_i = float(i);
        processor->reg[i].busy = 0;
    }

    // Initialize reservation stations based on configuration data
    int num_reservation_stations = configs.add_nr_reservation + configs.mul_nr_reservation + configs.div_nr_reservation;
    processor->reservation_stations = (reservation_station_t*)malloc(sizeof(reservation_station_t) * num_reservation_stations);
    if (processor->reservation_stations == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize reservation stations
    for (int i = 0; i < num_reservation_stations; i++)
    {
        init_reservation_station(&(processor->reservation_stations[i]), i + 1);
    }

    // Initialize units based on configuration data
    int num_units = configs.add_nr_units + configs.mul_nr_units + configs.div_nr_units;
    processor->units = (unit_t*)malloc(sizeof(unit_t) * num_units);
    if (processor->units == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize reservation stations
    uint32_t cdb_id = 0;
    for (int i = 0; i < num_units; i++)
    {
        if (i < configs.add_nr_units)
        {
            cdb_id = ADD_OPCODE;
        }
        else if (i < configs.add_nr_units + configs.mul_nr_units)
        {
            cdb_id = MUL_OPCODE;
        }
        else // i must be smaller than num_units, and up till there - the units are div units
        {
            cdb_id = DIV_OPCODE;
        }
        // TODO probably no need for the cdb id here
        init_unit(&(processor->units[i]), cdb_id);
    }

    // Initialize CDBs
    init_CDB(&(processor->cdb_add), ADD_OPCODE);
    init_CDB(&(processor->cdb_mul), MUL_OPCODE);
    init_CDB(&(processor->cdb_div), DIV_OPCODE);

    // TODO please explain to me the usage of memin and memin_size
    // Copy data from memin_data to memory locations in the processor

    // Additional initialization based on parsed data

    processor->instruction_que_size = 0;
    processor->instruction_que = (instruction_t*)malloc(sizeof(instruction_t)*INSTRUCTION_QUEUE_SIZE);
    if (processor->instruction_que == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }

    for (uint32_t i=0; i<INSTRUCTION_QUEUE_SIZE; ++i)
    {
        init_instruction(&(processor->instruction_que[i]));
    }

    processor->conf = configs;
    processor->halted = 0;
}

/**
 * @brief Prints the configurations.
 *
 * @param configurations The configurations to be printed.
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

/**
 * @brief Parses the configurations from a config file and updates the configurations structure.
 *
 * @param config_file Pointer to the config file.
 * @param configurations Pointer to the configurations structure to be updated.
 * @return 0 on success, -1 if the file is too big, 1 if memory allocation fails.
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

/**
 * @brief Parses the memory input file and updates the instruction array.
 *
 * @param memin_file Pointer to the memory input file.
 * @param instructions Pointer to the instruction array to be updated.
 * @return 0 on success, -1 if the file is too big, 1 if memory allocation fails.
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

        raw_instruction = strtoul(line, NULL, 16);
        instructions[instruction_index].pc = instruction_index;
        instructions[instruction_index].raw_instruction = raw_instruction;
        instructions[instruction_index].opcode = (raw_instruction & opcode_mask) >> 24;
        instructions[instruction_index].dst = (raw_instruction & dst_mask) >> 20;
        instructions[instruction_index].src0 = (raw_instruction & src0_mask) >> 16;
        instructions[instruction_index].src1 = (raw_instruction & src1_mask) >> 12;
        printf("opcode: %u; dst: %u; src0: %u; src1: %u\n", instructions[instruction_index].opcode, instructions[instruction_index].dst, instructions[instruction_index].src0, instructions[instruction_index].src1);
        ++instruction_index;

        memset(line, 0, sizeof(line));

        file_contents += content_index;
        if (file_contents_allocated_buffer + file_size <= file_contents)
        {
            break;
        }
    }

    free(file_contents_allocated_buffer);
    file_contents_allocated_buffer = NULL;
    return 0;
}


/**
 * @brief Fetches instructions from the instruction array and adds them to the instruction queue.
 *
 * @param instructions Pointer to the instruction array.
 * @param processor Pointer to the processor structure.
 */
void fetch_instructions(instruction_t** instructions, processor_t* processor)
{
    printf("que_size: %d\n", processor->instruction_que_size);

    // Instruction queue is full - nothing to do
    if (processor->instruction_que_size == INSTRUCTION_QUEUE_SIZE)
    {
        return;
    }
    // Fetching HALT instruction and station are not empty - skip
    if (processor->instruction_que_size != 0 && processor->instruction_que[processor->instruction_que_size - 1].opcode == HALT_OPCODE)
    {
        return;
    }

    // Only one place in the instruction queue - fetch only one instruction
    processor->instruction_que[processor->instruction_que_size] = *instructions[0];
    ++processor->instruction_que_size;
    // last instruction in memory
    if (instructions[0]->opcode == HALT_OPCODE)
    {
        return;
    }
    ++(*instructions);

    if (processor->instruction_que_size < INSTRUCTION_QUEUE_SIZE - 1)
    {
        // At least two places in the instruction queue - fetch two instructions
        processor->instruction_que[processor->instruction_que_size] = *instructions[0];
        ++processor->instruction_que_size;
        // last instruction in memory
        if (instructions[0]->opcode == HALT_OPCODE)
        {
            return;
        }
        ++(*instructions);
    }
    //TODO MAYBE_BUG why is PC not copied? is it it needed?
}

/**
 * @brief Executes the instructions based on the processor's configurations and reservation stations.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 */
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
        if (processor->reservation_stations[i].busy == 0 || 
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
            // TODO check that this works cause it doesn't seem that way
            if (
                !((opcode == ADD_OPCODE || opcode == SUB_OPCODE) && unit_id == ADD_OPCODE) &&
                !(opcode == MUL_OPCODE && unit_id == MUL_OPCODE) &&
                !(opcode == DIV_OPCODE && unit_id == DIV_OPCODE)
                )
            {
                continue;
            }
            
            processor->units[j].busy = 1;
            break;
        }
        
        processor->reservation_stations[i].start_cycle = cycle;
    }
}

/**
 * @brief Executes the task if it has finished and returns the result.
 *
 * @param cycle The current clock cycle.
 * @param reservation_station The reservation station containing the task to be executed.
 * @param conf The configurations of the processor.
 * @param halted Pointer to a flag indicating whether the processor is halted.
 * @param cdb_add Pointer to the ADD Common Data Bus (CDB).
 * @param cdb_mul Pointer to the MUL Common Data Bus (CDB).
 * @param cdb_div Pointer to the DIV Common Data Bus (CDB).
 * @return 1 if the task has finished, 0 otherwise.
 */
uint8_t execute_task_if_finished_and_get_result(uint32_t cycle, reservation_station_t reservation_station, configurations_t conf, uint8_t* halted, cdb_t* cdb_add, cdb_t* cdb_mul, cdb_t* cdb_div)
{
    long delay = 0;
    cdb_t* cdb = NULL;
    float result = 0;

    switch (reservation_station.op)
    {
    // TODO warning converting to float
    case ADD_OPCODE:
        cdb = cdb_add;
        result = reservation_station.vj + reservation_station.vk;
        delay = conf.add_delay;
        break;
    case SUB_OPCODE:
        cdb = cdb_add;
        result = reservation_station.vj - reservation_station.vk;
        delay = conf.add_delay;
        break;
    case MUL_OPCODE:
        cdb = cdb_mul;
        result = reservation_station.vj * reservation_station.vk;
        delay = conf.mul_delay;
        break;
    case DIV_OPCODE:
        if (reservation_station.vk == 0)
        {
            printf("Division by zero!\n");
            exit(-1);
        }
        result = reservation_station.vj / reservation_station.vk;
        cdb = cdb_div;
        delay = conf.div_delay;
        break;
    case HALT_OPCODE:
        // TODO this should never happen - make sure of it
        *halted = 1;
        printf("\n\nCheck that halted is working!\n\n");
        exit(13);
        break;
    default:
        exit(11);
        break;
    }

    // TODO check that >= is not buggy
    if (cycle >= reservation_station.start_cycle + delay && cdb->busy == 0)
    {
        cdb->busy = 1;
        cdb->pc = reservation_station.ins;
        cdb->output = result;
        cdb->station_id = reservation_station.station_id;
        return 1;
    }

    return 0;
}

/**
 * @brief Ends the execution of the processor if possible.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 */
void end_execution_if_possible(processor_t *processor, uint32_t cycle)
{
    // TODO Wait for the cdb to be free!
    uint32_t nr_reservation = 
        processor->conf.add_nr_reservation + 
        processor->conf.mul_nr_reservation + 
        processor->conf.div_nr_reservation;
    uint32_t nr_units =
        processor->conf.add_nr_units +
        processor->conf.mul_nr_units +
        processor->conf.div_nr_units;
        
    for (uint32_t i = 0; i < nr_reservation; ++i)
    {
        if (processor->reservation_stations[i].busy == 0 || processor->reservation_stations[i].start_cycle == 0)
        {
            continue;
        }

        if (execute_task_if_finished_and_get_result(cycle, processor->reservation_stations[i], processor->conf, &(processor->halted), &(processor->cdb_add), &(processor->cdb_mul), &(processor->cdb_div)) == 0)
        {
            continue;
        }
        
        // TODO write result and station id (or something else?) somewhere.
    }

    // TODO maybe we can return here false cause the res stations should never be empty
}

/**
 * @brief Writes the Common Data Bus (CDB) with the results of finished tasks.
 *
 * This function writes the results of finished tasks to the Common Data Bus (CDB) in the processor. 
 * It iterates through all reservation stations to check if any task has finished execution in the current cycle. 
 * If a task has finished, it retrieves the result and updates dependent reservation stations, registers, and execution units accordingly. 
 * Finally, it marks the execution unit as not busy, empties the reservation station, and prints debugging information.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 * @param cdb Pointer to the CDB structure to be updated.
 */
void write_cdb(processor_t* processor, uint32_t cycle, cdb_t* cdb,FILE* tracecdb_file)
{
    uint32_t nr_reservation =
        processor->conf.add_nr_reservation +
        processor->conf.mul_nr_reservation +
        processor->conf.div_nr_reservation;
    uint32_t nr_units =
        processor->conf.add_nr_units +
        processor->conf.mul_nr_units +
        processor->conf.div_nr_units;

    const float result = cdb->output;
    const uint32_t sending_station_id = cdb->station_id;
    const uint32_t pc = cdb->pc;

    if (sending_station_id == 0)
    {
        printf("Weird!\n\n");
        exit(14);
    }

    // Go over awaiting reservation stations and update with result
    for (uint32_t j = 0; j < nr_reservation; ++j)
    {
        if (processor->reservation_stations[j].qj == sending_station_id)
        {
            processor->reservation_stations[j].vj = result;
            processor->reservation_stations[j].qj = 0;
        }

        if (processor->reservation_stations[j].qk == sending_station_id)
        {
            processor->reservation_stations[j].vk = result;
            processor->reservation_stations[j].qk = 0;
        }
    }

    // Go over awaiting registers and update with result
    for (uint32_t k = 0; k < NUM_REGISTERS; k++)
    {
        if (processor->reg[k].q_i == sending_station_id)
        {
            processor->reg[k].v_i = result;
            processor->reg[k].q_i = 0;
        }
    }

    // mark execution unit as not busy
    for (uint32_t j = 0; j < nr_units; j++)
    {
        if (processor->units[j].pc == pc)
        {
            processor->units[j].busy = 0;
            break;
        }
    }

    // empty reservation station
    printf("insturction %d, set reservation station %d to NOT busy\n", pc, sending_station_id);

    // TODO should we do this here?
    init_reservation_station(&(processor->reservation_stations[sending_station_id-1]), sending_station_id);

    printf("write_cdb: cycle=%d, result=%f from station id %d from instruction number %d\n", cycle, result, sending_station_id, pc);
    write_tracecdb_file(tracecdb_file, cycle, result, sending_station_id, pc);
}

/**
 * @brief Writes the Common Data Bus (CDB) with the results of finished tasks if possible.
 *
 * This function writes the results of finished tasks to the Common Data Bus (CDB) in the processor if the CDB is not busy. 
 * It checks if the CDB is busy and writes the result if it is not. 
 * It also updates the cycle of the CDB and marks the CDB as busy.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 */
void write_cdb_if_possible(processor_t* processor, uint32_t cycle, FILE* tracecdb_file)
{
    if (processor->cdb_add.busy == 1)
    {
        write_cdb(processor, cycle, &(processor->cdb_add), tracecdb_file);
        processor->cdb_add.busy = 0;
        // TODO is this correct?
        processor->cdb_add.cycle = cycle;
    }

    if (processor->cdb_mul.busy == 1)
    {
        write_cdb(processor, cycle, &(processor->cdb_mul), tracecdb_file);
        processor->cdb_mul.busy = 0;
        processor->cdb_mul.cycle = cycle;
    }

    if (processor->cdb_div.busy == 1)
    {
        write_cdb(processor, cycle, &(processor->cdb_div), tracecdb_file);
        processor->cdb_div.busy = 0;
        processor->cdb_div.cycle = cycle;
    }

    //if (processor->cdb_mul.busy == 0 && result_mul_exists == 1)
    //{
    //    processor->cdb_mul.output = result_mul;
    //    processor->cdb_mul.busy = 1;
    //}
    //else if (processor->cdb_mul.busy == 1)
    //{
    //    processor->cdb_mul.busy = 0;
    //}

    //if (processor->cdb_div.busy == 0 && result_div_exists == 1)
    //{
    //    processor->cdb_div.output = result_div;
    //    processor->cdb_div.busy = 1;
    //}
    //else if (processor->cdb_div.busy == 1)
    //{
    //    processor->cdb_div.busy = 0;
    //}

}

/**
 * @brief Gets the start and end indices of the reservation stations based on the opcode.
 *
 * @param opcode The opcode of the instruction.
 * @param configurations The configurations of the processor.
 * @param o_start_index Pointer to store the start index of the reservation stations.
 * @param o_end_index Pointer to store the end index of the reservation stations.
 * @param halted Pointer to a flag indicating whether the processor is halted.
 */
void get_reservation_stations_indices_by_opcode(uint32_t opcode, configurations_t configurations, uint32_t *o_start_index, uint32_t *o_end_index, uint8_t* halted)
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
        *o_start_index = 0;
        *o_end_index = 0;
        *halted = 1;
        break;
    default:
        printf("try_to_issue_instructions: Got invalid instruction: %d\n", opcode);
        exit(-1);
    }
}

/**
 * @brief Issues instructions to the reservation stations based on the processor's instruction queue.
 *
 * @param processor Pointer to the processor structure.
 * @param o_first_instruction Pointer to the first instruction to be issued.
 * @param o_second_instruction Pointer to the second instruction to be issued.
 * @param current_cycle The current clock cycle.
 */
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
    if (processor->instruction_que_size == 0 || processor->instruction_que_size == 1 && processor->halted == 1)
    {
        return;
    }

    // TODO - refactor

    uint32_t start_point = 0;
    uint32_t end_point = processor->conf.add_nr_reservation;
    uint32_t num_instructions_issued = 0;

    get_reservation_stations_indices_by_opcode(processor->instruction_que[0].opcode, processor->conf, &start_point, &end_point, &(processor->halted));
    for (uint32_t i = start_point; i < end_point; ++i)
    {
        if (processor->reservation_stations[i].busy == 0)
        {
            *o_first_instruction = processor->instruction_que[0];
            printf("insturction %d, set reservation station %d to busy\n", o_first_instruction->pc, processor->reservation_stations[i].station_id);
            processor->reservation_stations[i].busy = 1;
            o_first_instruction->reservation_station = &processor->reservation_stations[i];
            ++num_instructions_issued;
            o_first_instruction->cycle_issued = current_cycle;
            break;
        }
    }

    // there was only one instruction to issue, don't issue another one
    if (processor->instruction_que_size != 1)
    {
        get_reservation_stations_indices_by_opcode(processor->instruction_que[1].opcode, processor->conf, &start_point, &end_point, &(processor->halted));
        for (uint32_t i = start_point; i < end_point; ++i)
        {
            if (processor->reservation_stations[i].busy == 0)
            {
                *o_second_instruction = processor->instruction_que[1];
                printf("insturction %d, set reservation station %d to busy\n", o_second_instruction->pc, processor->reservation_stations[i].station_id);
                processor->reservation_stations[i].busy = 1;
                o_second_instruction->reservation_station = &processor->reservation_stations[i];
                ++num_instructions_issued;
                o_second_instruction->cycle_issued = current_cycle;
                break;
            }
        }
    }

    memmove(processor->instruction_que, processor->instruction_que + num_instructions_issued, (INSTRUCTION_QUEUE_SIZE - num_instructions_issued)*sizeof(instruction_t));

    processor->instruction_que_size -= num_instructions_issued;

    // TODO maybe replace with for loop
    memset(processor->instruction_que + processor->instruction_que_size, 0, num_instructions_issued * sizeof(instruction_t));
}

/**
 * @brief Assigns register values to the reservation station if possible.
 *
 * This function assigns register values to the reservation station if the required registers are ready. 
 * It checks if the source registers are not busy and assigns their values to the reservation station. 
 * If the source registers are busy, it assigns the reservation station's tag to the source register's tag. Finally, it marks the destination register as busy.
 *
 * @param regs Array of registers to be updated.
 * @param instruction Instruction to be assigned to the reservation station.
 */
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


/**
 * @brief Assigns an instruction to a reservation station and updates its parameters.
 *
 * This function checks if the given instruction has been assigned to a reservation station. If not, it skips the assignment. Otherwise, it updates the reservation station's parameters such as the instruction's program counter and opcode. Additionally, it calls the helper function `assign_register_values_if_possible` to update the register values in the reservation station if the required registers are ready.
 *
 * @param processor Pointer to the processor structure.
 * @param instruction Instruction to be assigned to a reservation station.
 */
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

/**
 * @brief Sets pending registers to busy status.
 *
 * This function iterates through all registers and updates their busy status. If a register is currently not busy (busy = 0), it remains unchanged. If a register is already busy (busy = 1), its busy status is updated to indicate pending status (busy = 2).
 *
 * @param regs Array of registers to be updated.
 */
void set_pending_register_to_busy(register_t *regs)
{
    for (uint32_t i = 0; i < NUM_REGISTERS; ++i)
    {
        regs[i].busy = (regs[i].busy == 0) ? 0 : 2;
    }
}


/**
 * @brief Checks if the processor can exit the simulation.
 *
 * This function checks whether all reservation stations in the processor are empty, indicating that there are no pending instructions to execute. It iterates through all reservation stations and returns 1 if all stations are empty, indicating that the processor can exit the simulation. Otherwise, it returns 0.
 *
 * @param reservation_stations Array of reservation stations.
 * @param config Configuration settings for the processor.
 * @return 1 if the processor can exit, 0 otherwise.
 */
uint8_t check_if_can_exit(reservation_station_t* reservation_stations, configurations_t config)
{
    uint32_t nr_reservation = config.add_nr_reservation + config.mul_nr_reservation + config.div_nr_reservation;
    for (uint32_t i = 0; i < nr_reservation; i++)
    {
        // return 0 (false) if there is a non-empty reservation station
        if (reservation_stations[i].busy != 0)
        {
            return 0;
        }
    }

    // if we reached here, all reservation stations are empty
    return 1;
}



/**
 * @brief Runs the processor to execute instructions.
 *
 * This function runs the processor in a loop to execute instructions. 
 * It continuously fetches, issues, and executes instructions until the processor is halted and there are no more instructions to issue. 
 * Internally, it orchestrates the execution flow by calling various helper functions to perform different stages of instruction execution such as starting execution if possible, 
 * writing to the common data bus (CDB), fetching instructions, issuing instructions to reservation stations, and setting pending registers to busy status.
 *
 * @param processor Pointer to the processor structure.
 * @param instructions Array of instructions to be executed.
 */
void run_processor(processor_t* processor, instruction_t* instructions)
{
    //open the file to write the traceinst file
    FILE* tracecdb_file = fopen("tracecdb.txt", "w");
    if (tracecdb_file == NULL)
    {
        printf("Failed to open traceinst file!\n");
        exit(1);
    }

    uint32_t cycle = 0;

    instruction_t first_instruction, second_instruction;
    init_instruction(&first_instruction);
    init_instruction(&second_instruction);

    while (1)
    {
        // Internally, all these functions do their job only when everything is ready for it
        
        // Should we first fetch and then issue? can it happen the other way around?

        // The reservation stations list will be constructed like this:
            // first add_nr_res_stations are the adder stations
            // the mul_nr_res_stations are the multiplier stations
            // the the last div_nr_res_stations are the divider stations
        write_cdb_if_possible(processor, cycle, tracecdb_file);

        start_execution_if_possible(processor, cycle);
        
        end_execution_if_possible(processor, cycle);

        //get_instructions_to_execute(processor);
        //execute(processor);
        
        //write_tracecdb_file(cycle);
    
        // if got halted, continue execution loop but don't fetch and issue new instructions.
        // loop until there are no more instructions to issue
        if (processor->halted == 1)
        {
            ++cycle;
            if (check_if_can_exit(processor->reservation_stations, processor->conf) == 1)
            {
                break;
            }
            continue;
        }

        fetch_instructions(&instructions, processor);
        issue_instructions(processor, &first_instruction, &second_instruction, cycle);

        assign_instruction_to_reservation_station(processor, first_instruction);
        assign_instruction_to_reservation_station(processor, second_instruction);

        init_instruction(&first_instruction);
        init_instruction(&second_instruction);

        set_pending_register_to_busy(processor->reg);
        ++cycle;
    }
    //open file tranceinst to write the traceinst file
    FILE* traceinst_file = fopen("traceinst.txt", "w");
    if (traceinst_file == NULL)
    {
        printf("Failed to open traceinst file!\n");
        exit(1);
    }
    //open file regout to write the regout file
    FILE* regout_file = fopen("regout.txt", "w");
    if (regout_file == NULL)
    {
        printf("Failed to open regout file!\n");
        exit(1);
    }
    //write the traceinst file
    write_traceinst_file(traceinst_file, instructions);
    //write the regout file
    write_regout_file(regout_file, processor->reg);
    //close the files
    fclose(traceinst_file);
    fclose(regout_file);
    fclose(tracecdb_file);
}

int write_regout_file(FILE* regout_file, register_t* reg) {
    // Implementation of writing the register output file
    for (size_t i = 0; i < 16; i++)
    {
        fprintf(regout_file, "%f\n", reg[i].v_i);
    }
    
    // Return 0 on success, non-zero value on failure
    return 0;
}


// do we need it?
char *find_tag_name(char *tag, instruction_t instruction)
{
    if (instruction.opcode == ADD_OPCODE || instruction.opcode == SUB_OPCODE)
    {
        sprintf(tag, "ADD%d", instruction.reservation_station->station_id);
    }
    else if (instruction.opcode == MUL_OPCODE)
    {
        sprintf(tag, "MUL%d", instruction.reservation_station->station_id);
    }
    else if (instruction.opcode == DIV_OPCODE)
    {
        sprintf(tag, "DIV%d", instruction.reservation_station->station_id);
    }
    else
    {
        sprintf(tag, "HALT");
    }
    return tag;
}


int write_traceinst_file(FILE* traceinst_file, instruction_t* instructions) {
    // Implementation of writing the trace instruction file
    int i = 0;
    char tag_name[10] = {}; // Add initializer
    //find the tag name in char while getting intruction.tag[i].tag

    while (instructions[i].raw_instruction != 0 )
    {
        if (instructions[i].opcode == HALT_OPCODE)
        {
            return 0;
        }
        else
        {
            find_tag_name(tag_name, instructions[i]);
            fprintf(traceinst_file, "%08x %d %s %d %d %d %d\n", instructions[i].raw_instruction,instructions[i].pc,tag_name, instructions[i].cycle_issued, instructions[i].cycle_execute_start, instructions[i].cycle_execute_end,instructions[i].cycle_cdb);
        }
        i++;
    }
    // Return 0 on success, non-zero value on failure
    return 0;
}

// need to run each cycle when using CDB
int write_tracecdb_file(FILE* tracecdb_file, uint32_t cycle,uint32_t result,uint32_t sending_station_id,uint32_t pc) {
    // Implementation of writing the trace CDB file
    fprintf(tracecdb_file, "%d %d %f %d\n", cycle, pc, result, sending_station_id);
    // Return 0 on success, non-zero value on failure
    return 0;
}

/**
 * @brief Main function to run the processor simulation.
 *
 * This main function is responsible for initializing the processor, parsing configuration and input files, and executing the processor simulation by calling the `run_processor` function. It takes command-line arguments for configuration file and input file paths and performs error handling for file opening and parsing. After running the processor simulation, it returns 0 to indicate successful execution.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return 0 if execution is successful, 1 otherwise.
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
    configurations_t configs;

    errno_t err_config_file = 0;
    errno_t err_memin_file = 0;
    errno_t err_regout_file = 0;
    errno_t err_traceinst_file = 0;
    errno_t err_tracecdb_file = 0;

    processor_t processor;

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

    init_processor(&processor, configs);
    run_processor(&processor, instructions);
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
