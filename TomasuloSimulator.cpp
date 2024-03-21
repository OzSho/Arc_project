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
    configurations->add_nr_units = 0;
    configurations->mul_nr_units = 0;
    configurations->div_nr_units = 0;
    configurations->add_nr_reservation = 0;
    configurations->mul_nr_reservation = 0;
    configurations->div_nr_reservation = 0;
    configurations->add_delay = 0;
    configurations->mul_delay = 0;
    configurations->div_delay = 0;
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
void init_processor(processor_t* processor, configurations_t configs)
{
    // Initialize processor state variables
    
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
    for (int i = 0; i < num_units; i++)
    {
        // TODO must resolve cdb id
        init_unit(&(processor->units[i]), 1);
    }

    // Initialize CDBs
    init_CDB(&(processor->cdb_add), 1);
    init_CDB(&(processor->cdb_mul), 2);
    init_CDB(&(processor->cdb_div), 3);

    // TODO please explain to me the usage of memin and memin_size
    // Copy data from memin_data to memory locations in the processor

    // Additional initialization based on parsed data

    processor->instruction_que_size = 0;
    processor->instruction_que = (instruction_t*)malloc(sizeof(instruction_t)*INSTRUCTION_QUEUE_SIZE);
    for (uint32_t i=0; i<INSTRUCTION_QUEUE_SIZE; ++i)
    {
        init_instruction(&(processor->instruction_que[i]));
    }

    processor->conf = configs;
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

int write_regout_file() {
    // Implementation of writing the register output file
    // Return 0 on success, non-zero value on failure
    return 0;
}

int write_traceinst_file() {
    // Implementation of writing the trace instruction file
    // Return 0 on success, non-zero value on failure
    return 0;
}

int write_tracecdb_file() {
    // Implementation of writing the trace CDB file
    // Return 0 on success, non-zero value on failure
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

    // Check if the instruction queue is full
    if (processor->instruction_que_size == INSTRUCTION_QUEUE_SIZE)
    {
        return; // Nothing to do, instruction queue is full
    }

    // Fetch the first instruction and add it to the instruction queue
    processor->instruction_que[processor->instruction_que_size] = *instructions[0];
    ++processor->instruction_que_size;
    ++(*instructions);

    // Check if there is enough space in the instruction queue for another instruction
    if (processor->instruction_que_size < INSTRUCTION_QUEUE_SIZE - 1)
    {
        // Fetch the second instruction and add it to the instruction queue
        processor->instruction_que[processor->instruction_que_size] = *instructions[0];
        ++processor->instruction_que_size;
        ++(*instructions);
    }
    // TODO: Consider copying the PC (Program Counter) to the instruction queue
}

/**
 * @brief Executes the instructions based on the processor's configurations and reservation stations.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 */
void start_execution_if_possible(processor_t* processor, uint32_t cycle) {
    // Implementation of executing the instructions based on configurations and reservation stations
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
            if (
                !((opcode == ADD_OPCODE || opcode == SUB_OPCODE) && unit_id == ADD_ID) &&
                !(opcode == MUL_OPCODE && unit_id == MUL_ID) &&
                !(opcode == DIV_OPCODE && unit_id == DIV_ID)
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
 * @param result Pointer to store the result of the task.
 * @param halted Pointer to a flag indicating whether the processor is halted.
 * @return 1 if the task has finished, 0 otherwise.
 */
uint8_t execute_task_if_finished_and_get_result(uint32_t cycle, reservation_station_t reservation_station, configurations_t conf, float* result, uint8_t* halted)
{
    long delay = 0;
    
    switch (reservation_station.op)
    {
    // TODO warning converting to float
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
        *halted = 1;
        printf("\n\nCheck that halted is working!\n\n");
    default:
        exit(11);
        break;
    }

    return (cycle == reservation_station.start_cycle + delay) ? 1 : 0;
}

/**
 * @brief Writes the Common Data Bus (CDB) with the results of finished tasks.
 *
 * This function writes the results of finished tasks to the Common Data Bus (CDB) in the processor. It iterates through all reservation stations to check if any task has finished execution in the current cycle. If a task has finished, it retrieves the result and updates dependent reservation stations, registers, and execution units accordingly. Finally, it marks the execution unit as not busy, empties the reservation station, and prints debugging information.
 *
 * @param processor Pointer to the processor structure.
 * @param cycle The current clock cycle.
 */
void write_cdb(processor_t *processor, uint32_t cycle) {
    // Calculate the total number of reservation stations and units
    uint32_t nr_reservation = processor->conf.add_nr_reservation +
                              processor->conf.mul_nr_reservation +
                              processor->conf.div_nr_reservation;
    uint32_t nr_units = processor->conf.add_nr_units +
                        processor->conf.mul_nr_units +
                        processor->conf.div_nr_units;

    float result = 0;

    // Iterate through all reservation stations
    for (uint32_t i = 0; i < nr_reservation; ++i) {
        // Skip if the station is not busy or hasn't started execution
        if (processor->reservation_stations[i].busy == 0 || processor->reservation_stations[i].start_cycle == 0) {
            continue;
        }

        // Check if task has finished execution and get the result
        if (execute_task_if_finished_and_get_result(cycle, processor->reservation_stations[i], processor->conf, &result, &(processor->halted)) == 0) {
            continue; // Skip if task hasn't finished
        }

        // Update dependent reservation stations with the result
        for (uint32_t j = 0; j < nr_reservation; ++j) {
            if (i == j) {
                continue; // Skip the current station
            }

            if (processor->reservation_stations[j].qj == processor->reservation_stations[i].station_id) {
                processor->reservation_stations[j].vj = result; // Update Vj if waiting for the result
            }

            if (processor->reservation_stations[j].qk == processor->reservation_stations[i].station_id) {
                processor->reservation_stations[j].vk = result; // Update Vk if waiting for the result
            }
        }

        // Update awaiting registers with the result
        for (uint32_t k = 0; k < NUM_REGISTERS; k++) {
            if (processor->reg[k].q_i == processor->reservation_stations[i].station_id) {
                processor->reg[k].v_i = result; // Update Vi if waiting for the result
            }
        }

        // Store cycle_cdb in the instruction struct (TODO: where to store this?)
        

        instruction_t instruction;
        instruction.cycle_cdb = cycle;

        // Mark the execution unit as not busy
        for (uint32_t j = 0; j < nr_units; j++) {
            if (processor->reservation_stations[i].ins == processor->units[j].pc) {
                processor->units[j].busy = 0; // Mark the unit as not busy
                break;
            }
        }

        // Empty the reservation station
        processor->reservation_stations[i].busy = 0;
        init_reservation_station(&(processor->reservation_stations[i]), processor->reservation_stations[i].station_id); // Reinitialize the station

        // Print debugging information
        printf("write_cdb: cycle=%d, result=%d\n", cycle, result);
    }
}


/**
 * @brief Gets the start and end indices of the reservation stations based on the opcode.
 *
 * @param opcode The opcode of the instruction.
 * @param configurations The configurations of the processor.
 * @param o_start_index Pointer to store the start index of the reservation stations.
 * @param o_end_index Pointer to store the end index of the reservation stations.
 */
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
        // TODO: what should we do here? probably nothing
        break;
    default:
        printf("try_to_issue_instructions: Got invalid instruction: %d\n", opcode);
        exit(-1);
    }
}

/**
 * @brief Issues 1 instruction from the instruction queue to the reservation stations.
 *
 * This function attempts to issue at most one instruction from the instruction queue to the reservation stations in the processor. It marks the newly occupied reservation stations as busy and writes the station id of the reservation station inside the instruction.
 * After issuing instructions, it reorganizes the instruction queue by removing the issued instructions.
 *
 * @param processor Pointer to the processor structure.
 * @param o_instruction Pointer to store the issued instruction.
 * @param index Index of the instruction in the instruction queue.
 * @param current_cycle Current cycle of the processor.
 * @return 1 if an instruction was issued, 0 otherwise.
 */
int issue_instruction(processor_t* processor, instruction_t* o_instruction, uint32_t index, const uint32_t current_cycle) {
    if (processor == NULL) {
        printf("issue_instruction: Got NULL processor!\n");
        return 0;
    }

    if (processor->instruction_que_size <= index) {
        return 0;
    }

    uint32_t start_point = 0;
    uint32_t end_point = processor->conf.add_nr_reservation;

    get_reservation_stations_indices_by_opcode(processor->instruction_que[index].opcode, processor->conf, &start_point, &end_point);

    for (uint32_t i = start_point; i < end_point; ++i) {
        if (processor->reservation_stations[i].busy == 0) {
            *o_instruction = processor->instruction_que[index];
            processor->reservation_stations[i].busy = 1;
            o_instruction->reservation_station = &processor->reservation_stations[i];
            o_instruction->cycle_issued = current_cycle;
            return 1;
        }
    }

    return 0;
}

/**
 * @brief Issues instructions from the instruction queue to the reservation stations.
 *
 * This function attempts to issue at most two instructions from the instruction queue to the reservation stations in the processor. It marks the newly occupied reservation stations as busy and writes the station id of the reservation station inside the instruction.
 * After issuing instructions, it reorganizes the instruction queue by removing the issued instructions.
 *
 * @param processor Pointer to the processor structure.
 * @param o_first_instruction Pointer to store the first issued instruction.
 * @param o_second_instruction Pointer to store the second issued instruction.
 * @param current_cycle Current cycle of the processor.
 */
void issue_instructions(processor_t* processor, instruction_t* o_first_instruction, instruction_t* o_second_instruction, const uint32_t current_cycle) {
    if (processor == NULL) {
        printf("issue_instructions: Got NULL processor!\n");
        return;
    }

    if (processor->instruction_que_size == 0) {
        return;
    }

    uint32_t num_instructions_issued = 0;

    if (issue_instruction(processor, o_first_instruction, 0, current_cycle)) {
        num_instructions_issued++;
    }

    if (processor->instruction_que_size > 1 && issue_instruction(processor, o_second_instruction, 1, current_cycle)) {
        num_instructions_issued++;
    }

    memmove(processor->instruction_que, processor->instruction_que + num_instructions_issued, (INSTRUCTION_QUEUE_SIZE - num_instructions_issued)*sizeof(instruction_t));
    processor->instruction_que_size -= num_instructions_issued;
    memset(processor->instruction_que + processor->instruction_que_size, 0, num_instructions_issued * sizeof(instruction_t));
}

/**
 * @brief Updates the register values in the reservation station if possible.
 *
 * This function checks the status of source registers and updates the corresponding values and queue flags in the reservation station based on their availability. If a source register is not busy, its value is directly assigned to the reservation station; otherwise, the queue flag is updated with the register's queue id, and the value is set to UINT16_MAX to indicate that it's not available.
 *
 * @param regs Array of registers.
 * @param instruction Instruction containing source and destination register indices.
 */
void assign_register_values_if_possible(register_t* regs, instruction_t instruction) {
    update_register_value(regs, instruction.src0, &instruction.reservation_station->vj, &instruction.reservation_station->qj);
    update_register_value(regs, instruction.src1, &instruction.reservation_station->vk, &instruction.reservation_station->qk);

    regs[instruction.dst].busy = 1;
    regs[instruction.dst].q_i = instruction.reservation_station->station_id;
}

/**
 * @brief Updates the register value and queue flag in the reservation station based on the source register's status.
 *
 * This function updates the reservation station's value and queue flag based on the status of the specified source register.
 *
 * @param regs Array of registers.
 * @param src_index Index of the source register.
 * @param value Pointer to the value in the reservation station.
 * @param queue_flag Pointer to the queue flag in the reservation station.
 */
void update_register_value(register_t* regs, uint32_t src_index, uint32_t* value, uint32_t* queue_flag) {
    if (regs[src_index].busy != 2) {
        *value = regs[src_index].v_i;
        *queue_flag = 0;
    } else {
        *queue_flag = regs[src_index].q_i;
        *value = UINT16_MAX;
    }
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
void set_pending_register_to_busy(register_t *regs) {
    for (uint32_t i = 0; i < NUM_REGISTERS; ++i) {
        regs[i].busy = (regs[i].busy == 0) ? 0 : 2; // Set pending registers to busy
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
uint8_t check_if_can_exit(reservation_station_t* reservation_stations, configurations_t config) {
    uint32_t nr_reservation = config.add_nr_reservation + config.mul_nr_reservation + config.div_nr_reservation;
    for (uint32_t i = 0; i < nr_reservation; i++) {
        // Return 0 (false) if there is a non-empty reservation station
        if (reservation_stations[i].busy != 0) {
            return 0;
        }
    }

    // If we reached here, all reservation stations are empty
    return 1; // Return 1 (true) indicating the processor can exit
}


/**
 * @brief Runs the processor to execute instructions.
 *
 * This function runs the processor in a loop to execute instructions. It continuously fetches, issues, and executes instructions until the processor is halted and there are no more instructions to issue. Internally, it orchestrates the execution flow by calling various helper functions to perform different stages of instruction execution such as starting execution if possible, writing to the common data bus (CDB), fetching instructions, issuing instructions to reservation stations, and setting pending registers to busy status.
 *
 * @param processor Pointer to the processor structure.
 * @param instructions Array of instructions to be executed.
 */
void run_processor(processor_t* processor, instruction_t* instructions)
{
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

        start_execution_if_possible(processor, cycle);

        write_cdb(processor, cycle);

        //get_instructions_to_execute(processor);
        //execute(processor);
        
        //write_tracecdb_file(cycle);
    
        // if got halted, continue execution loop but don't fetch and issue new instructions.
        // loop until there are no more instructions to issue
        if (processor->halted == 1)
        {
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

    //write_traceinst_file(/*instructions*/);
    //write_regout_file(/*processor*/);
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
    // FILE *regout_file = NULL;
    // FILE *traceinst_file = NULL;
    // FILE *tracecdb_file = NULL;
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