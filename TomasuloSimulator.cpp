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
#define MAX_MEMIN_FILE_SIZE 40960
#define MAX_FILE_LINE_LEN 50
#define NUM_REGISTERS 16
#define ADD_ID 0
#define MUL_ID 1
#define DIV_ID 2

typedef struct {
    uint32_t cdb_id;        // Identifier for the CDB (ADD, MUL, DIV)
    uint32_t output;        // output transmitted on the CDB
    uint32_t busy;          // Flag indicating if the station is busy (1) or available (0)
} cdb_t;

typedef struct {
    uint32_t cdb_id;        // Identifier for the CDB (e.g., ADD, MUL, DIV)
    uint32_t tag;           // number of the CDB in the group
    uint32_t busy;          // Flag indicating if the station is busy (1) or available (0)
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
} reservation_station_t;

typedef struct {
    uint32_t opcode;            // Opcode of the instruction
    uint32_t dst;               // Destination register index
    uint32_t src0;              // Source register 0 index
    uint32_t src1;              // Source register 1 index
    uint32_t cycle_issued;      // Cycle when the instruction was issued
    uint32_t cycle_execute_start; // Cycle when execution of the instruction started
    uint32_t cycle_execute_end; // Cycle when execution of the instruction ended
    uint32_t cycle_cdb;         // Cycle when the result was written to the Common Data 
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
    reservation_station_t* reservation_stations; // Dynamically allocated array of reservation stations
    cdb_t cdb_add; // add Data Buses
    cdb_t cdb_mul; // mul Data Buses
    cdb_t cdb_div; // div Data Buses
    uint32_t cycle;
    configurations_t conf;
    unit_t* units; //  Dynamically allocated array of units
    instruction_t* instruction_que; // Array to store up to 16 instructions 
} processor_t;

//why do we need it?
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

void init_CDB(cdb_t* cdb,uint32_t cdb_id,uint32_t tag )        // output transmitted on the CDB
{
    //
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

    char line[MAX_FILE_LINE_LEN] = {0};

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
int parse_memin(FILE* memin_file)
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

    long command = 0;
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

        command = strtol(line, NULL, 16);
        opcode = (command & opcode_mask) >> 24;
        dst = (command & dst_mask) >> 20;
        src0 = (command & src0_mask) >> 16;
        src1 = (command & src1_mask) >> 12;
        printf("opcode: %u; dst: %u; src0: %u; src1: %u\n", opcode, dst, src0, src1);

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

int write_regout_file() {
    // Implementation of writing the register output file
    // Return 0 on success, non-zero value on failure
    return 0;
}

int write_traceinst_file() {
    // Implementation of writing the trace instruction file
    // Return 0 on success, non-zero value on failure
}

int write_tracecdb_file() {
    // Implementation of writing the trace CDB file
    // Return 0 on success, non-zero value on failure
}

void init_processor(processor_t* processor, configurations_t* configs, uint32_t* memin_data, size_t memin_size) {
    // Initialize processor state variables
    processor->pc = 0;
    
    // Initialize registers based on configuration data
    processor->reg = (float*)malloc(sizeof(float) * NUM_REGISTERS);
    if (processor->reg == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize registers to zero or other default values
    for(int i=0;i<16;i++)
    {
        processor->reg[i]=float(i);
    }
    // Initialize reservation stations based on configuration data
    int num_reservation_stations = configs->add_nr_reservation + configs->mul_nr_reservation + configs->div_nr_reservation;
    processor->reservation_stations = (reservation_station_t*)malloc(sizeof(reservation_station_t) * num_reservation_stations);
    if (processor->reservation_stations == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize reservation stations
    for(int i=0;i<num_reservation_stations;i++)
    {
        //
    }
    // Initialize units based on configuration data
    int num_units = configs->add_nr_units + configs->mul_nr_units + configs->div_nr_units;
    processor->units = (unit_t*)malloc(sizeof(unit_t) * num_units);
    if (processor->units == NULL) {
        // Handle memory allocation error
        exit(EXIT_FAILURE);
    }
    // Initialize CDBs
    for(int i=0;i<num_reservation_stations;i++)
    {
        //
    }
    // Copy data from memin_data to memory locations in the processor
    
    // Additional initialization based on parsed data
    
    // Free memory allocated for configurations after initialization
    free(configs);
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

    FILE *config_file = NULL;
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

    if (parse_memin(memin_file) != 0)
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
