#include <iostream>
#include <string>
#include <stdint.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <bitset>
#include <cmath>
#include <iomanip>
#include <vector>
#include <unordered_map>


#define MAX_PRIORITY 99999
#define INFINITE 999999999
typedef std::unordered_map<uint32_t, std::vector<int>> access_cycle_tracker_t ;


bool debug_mode = false;
bool output_to_file = false;
int max_debug_line = 1000000;

bool simulation = true;
bool graph = false;


// Enum for Inclusion Property
enum InclusionProperty {
    NON_INCLUSIVE = 0,
    INCLUSIVE = 1
};

// Enum for Replacement Policy
enum ReplacementPolicy {
    LRU = 0,
    FIFO = 1,
    OPTIMAL = 2
};

// Enum for Inclusion Property
enum memory_State {
    VALID = 0,
    INVALID = 1,
    DIRTY = 2
};


uint32_t hexAddrToUint32(std::string hexString);
std::string uint32ToHexChar(uint32_t value);


void printBinary(uint32_t value) {
    std::bitset<32> binary(value);  // Convert value to a 32-bit binary representation
    std::string binaryStr = binary.to_string();

    // Insert spaces every 4 bits
    for (int i = 0; i < 32; i += 4) {
        std::cout << binaryStr.substr(i, 4) << " ";
    }
    std::cout << std::endl;
}

uint32_t hexAddrToUint32(std::string hexString) {
    
    // Buffer to store the padded string if needed
    char paddedHex[9] = {0}; // 8 characters + 1 for null terminator

    // Determine the length of the input
    // int len = hexString.length();
    int zeros = 8 - hexString.length();
    zeros =  (zeros>=0) ? zeros : 0 ;

    int i = 0;
    while(i<zeros) paddedHex[i++] = '0';

    int j = 0;
    while(i <= 8) paddedHex[i++] = hexString[j++];    // number of zeros starts count from 1, but index from 0. So having n zeros means you need to start writing from nth index


    uint32_t address = std::stoul(paddedHex, nullptr, 16);
    return address;

}

std::string uint32ToHexChar_wzero(uint32_t value) {
    // Allocate space for 8 hex characters + null terminator
    char hexString[9];
    std::snprintf(hexString, sizeof(hexString), "%08x", value); // Format as zero-padded 8-digit hex
    return std::string(hexString);
}

std::string uint32ToHexChar(uint32_t value) {
    std::stringstream ss;
    ss << std::hex << std::nouppercase << value; // Convert to lowercase hex without padding
    return ss.str();
}


uint32_t get_address(uint32_t tag, uint32_t index, uint32_t offset, int n_index_bits, int n_offset_bits) {
    // Shift the tag left by (n_index_bits + offset_bits) and add the shifted index and offset.
    return (tag << (n_index_bits + n_offset_bits)) | (index << n_offset_bits) | offset;
}

uint32_t get_tag(uint32_t address, int n_index_bits, int n_offset_bits) {
    // Shift the address right by the number of index and offset bits to get the tag
    return address >> (n_index_bits + n_offset_bits);
}

uint32_t get_index(uint32_t address, int n_index_bits, int n_offset_bits) {
    // Shift right by the offset bits to discard them, then mask the index bits
    return (address >> n_offset_bits) & ((1 << n_index_bits) - 1);
}

uint32_t get_offset(uint32_t address, int n_offset_bits) {
    // Mask the offset bits
    return address & ((1 << n_offset_bits) - 1);
}


// Helper function to get policy name as a string
std::string get_policy_name(ReplacementPolicy policy) {
    switch (policy) {
        case LRU: return "LRU";
        case FIFO: return "FIFO";
        case OPTIMAL: return "OPTIMAL";
        default: return "Unknown";
    }
}

// Helper function to get inclusion name as a string
std::string get_inclusion_name(InclusionProperty inclusion) {
    switch (inclusion) {
        case INCLUSIVE: return "inclusive";
        case NON_INCLUSIVE: return "non-inclusive";
        default: return "unknown";
    }
}


access_cycle_tracker_t access_cycle_tracker ;

class CACHE {
private:
    std::string cache_level;
    int size;                 // Cache size in bytes
    int assoc;                // Associativity              // #columns in the matrix
    int block_size;           // Block size in bytes
    int set;                  // size/(assoc*block_size)    // #rows in the matrix
    ReplacementPolicy replacement_policy; // Replacement policy
    InclusionProperty inclusion_property; // Inclusion property

    int tagBits;
    int indexBits;
    int offsetBits;

    int n_reads, n_writes;
    int n_read_misses, n_write_misses;
    int n_write_backs, n_write_backs_to_main;
    double miss_rate;

    uint32_t * memory ;
    int* state_matrix ;
    int* priority_matrix ;
    int current_priority;   // aka current_cycle
    // int row = index/assoc;
    // int column = index%assoc;

    CACHE* lower_memory;
    CACHE* higher_memory;


public:
    // Constructor for L1 and optional L2 cache
    CACHE(std::string level, int size, int assoc, int block_size, ReplacementPolicy policy, InclusionProperty inclusion)
        : size(size), assoc(assoc), block_size(block_size), replacement_policy(policy), inclusion_property(inclusion), 
          n_reads(0), n_writes(0), n_read_misses(0), n_write_misses(0), n_write_backs(0), miss_rate(0), n_write_backs_to_main(0),
          lower_memory(NULL), higher_memory(NULL), current_priority(0), cache_level(level) {

        if(assoc == 0) set = 0;    
        else set = size / (assoc * block_size);

        offsetBits = std::log2(block_size);
        indexBits = std::log2(set);
        tagBits = 32 - indexBits - offsetBits;    // addressSize = 32

        memory = (uint32_t*) malloc (size*sizeof(uint32_t));
        state_matrix = (int*) malloc (size*sizeof(int));
        priority_matrix = (int*) malloc (size*sizeof(int));

        
        // std::cout << "Size is = " << size << "\n";
        // initialize the arrays
        for (int i=0; i<size; i++) {
            memory[i] = 0;
            state_matrix[i] = INVALID;
            priority_matrix[i] = -1;
        }
    }

    void free_memories() {
        free(memory);
        free(state_matrix);
        free(priority_matrix);
    }
    
    int get_size() {
        return size;
    }

    int get_write_backs_to_main() {
        return n_write_backs_to_main;
    }

    // Function to set higher and lower memory
    void set_memory_hierarchy(CACHE *higher, CACHE *lower) {
        lower_memory = lower;
        higher_memory = higher;
    }


    // Function to update cache state 
    void update_state(int index, memory_State state) {      // index is position in the memory array
        state_matrix[index] = state;
        // std::cout << "Updating cache state..." << std::endl;
    }

    
    int search_memory(uint32_t tag_bits, uint32_t index_bits) {
        int starting_index = index_bits*assoc;
        for(int i=0; i<assoc; i++) {
            if(memory[starting_index+i] == tag_bits) return (starting_index + i);
        }

        return -1; // address not found. MISS
    }

    int evict(uint32_t address, bool from_lower) {      // 0:successful;    1:error:address already removed;
        // if(!size || !set) return 2; // error code = 2 : cache doesnot exist

        uint32_t tag_bits = get_tag(address, indexBits, offsetBits); 
        uint32_t index_bits = get_index(address, indexBits, offsetBits);

        int index = search_memory(tag_bits, index_bits);
        if (index == -1) return 1;  // address already removed
        // L1 invalidated: 7b034500 (tag 3d81a2, index 16, dirty)
        if(debug_mode) {
            std::cout << cache_level << " invalidated: " << uint32ToHexChar(address)
                    << " (tag " << uint32ToHexChar(tag_bits) << ", index " << index_bits ; 
            if(state_matrix[index] == DIRTY) std::cout << ", dirty)\n";
            else std::cout << ", clean)\n";
        }


        if(state_matrix[index] == DIRTY) {      // write to main memory directly
            // L1 writeback to main memory directly
            n_write_backs_to_main += 1;     // n_write_backs += 1;
            if(debug_mode) {
                std::cout << cache_level << " writeback to main memory directly\n";
            }
        }    

        memory[index] = 0;
        state_matrix[index] = INVALID;
        priority_matrix[index] = -1;

        return 0;
    }


    // Function to allocate a new block 
    int allocate_block(uint32_t tag_bits, uint32_t index_bits) {    // tag_bits not needed here
        // if(!size || !set) return -1; // error code = -1 : cache doesnot exist
        // get_block_according_to_policy

        int set_start_block_index = index_bits*assoc;
        int victim_indx = set_start_block_index;

        if(replacement_policy == OPTIMAL) { 
            int current_cycle = current_priority;
            int max_diff = 0;

            for(int i=0; i<assoc; i++) {    // traverse through each assoc memory in an index

                if(state_matrix[set_start_block_index + i] == INVALID && memory[set_start_block_index + i] == 0) {    // empty now : memory[set_start_block_index + i] = 0;
                    victim_indx = set_start_block_index + i;
                    if(debug_mode) std::cout << cache_level << " victim: none\n" ;     // as state = INVALID and thus memory[index]==0
                    return victim_indx;
                    
                }
                
                // get addr and vector corresponding to that address. 
                uint32_t addr = get_address(memory[set_start_block_index + i], index_bits, 0, indexBits, offsetBits);

                std::vector<int>& access_cycle_vector = access_cycle_tracker[addr];

                // now get the lowest value which is greater than current cycle. (won't be equal, cause in that case it would not be allocating)
                int low_val = INFINITE;

                for (auto it = access_cycle_vector.begin(); it != access_cycle_vector.end(); ) {

                    if (*it > current_cycle && *it < low_val) low_val = *it;
                    
                    // // std::cout << "stored value = " << *it << "\t";
                    // int diff = *it - current_cycle;
                    // if (diff > max_diff) max_diff = diff;

                    // Remove elements that are less than or equal to the current cycle
                    if (*it <= current_cycle) {
                        it = access_cycle_vector.erase(it); // Erase and update iterator
                    } else {
                        ++it; // Move to next element if not erased
                    }
                    
                }

                int vect_diff = low_val - current_cycle;

                if(vect_diff > max_diff) {
                    max_diff = vect_diff;
                    victim_indx = set_start_block_index + i ;
                }

                // if(max_diff != 0 && max_diff >= priority_diff) {
                //     victim_indx = set_start_block_index + i;
                //     priority_diff = max_diff;

                // }
                    
            }
            

        } else {    // LRU & FIFO
            for(int i=1; i<assoc; i++) {  // set_start_block_index(i=0) already checked as lowest index
                if(priority_matrix[set_start_block_index + i] < priority_matrix[victim_indx]) 
                    victim_indx = set_start_block_index + i;
            }

        }

        // for debug print 
        if(debug_mode) {
            if(memory[victim_indx] != 0) {
                std::cout << cache_level << " victim: " << uint32ToHexChar(get_address(memory[victim_indx], index_bits, 0, indexBits, offsetBits))
                        << " (tag " << uint32ToHexChar(memory[victim_indx]) << ", index " << index_bits ; 
                if(state_matrix[victim_indx] == DIRTY) std::cout << ", dirty)\n";
                else std::cout << ", clean)\n";
            } else {    
                std::cout << cache_level << " victim: none\n" ;
            } 
        }


        uint32_t address = get_address(memory[victim_indx], index_bits, 0, indexBits, offsetBits);

        // take action according to memory state and inclusion property
        if(state_matrix[victim_indx] == DIRTY) {
            if(lower_memory)    // Dirty block :: write_back to lower
                lower_memory->write_back(address);  // writing the tag bits should not be enough. add(mask) set(index bits)

            n_write_backs += 1;
            state_matrix[victim_indx] = VALID;  // after writing back to lower memory, it should be valid.
        }

        if(inclusion_property == INCLUSIVE && higher_memory && memory[victim_indx] != 0) {   // remove from higher memory entity
            higher_memory->evict(address, true);    // writing the tag bits should not be enough. add(mask) set(index bits)
        }
        

        // allocate :: reset everything for the memory
        memory[victim_indx] = 0;
        priority_matrix[victim_indx] = -1;   // setting to lowest index. so that it is accessed. as not used yet
        state_matrix[victim_indx] = INVALID;    // it is usable now

        // std::cout << "returning\t";
        return victim_indx;

    }

    // Function to find and read a block
    int read(uint32_t address) {
        // if(!size || !set) return -1; // error code = -1 : cache doesnot exist

        if(debug_mode) std::cout << cache_level << " read : " << uint32ToHexChar(address) ; 
        n_reads += 1;

        uint32_t tag_bits = get_tag(address, indexBits, offsetBits);
        uint32_t index_bits = get_index(address, indexBits, offsetBits);

        if(debug_mode) std::cout << " (tag " << uint32ToHexChar(tag_bits) << ", index " << index_bits << ")\n" ; 

        // search the address 
        int index = search_memory(tag_bits, index_bits); 

        if(index == -1) {   // MISS. fetch from lower
            if(debug_mode) std::cout << cache_level << " miss\n";
            n_read_misses += 1;
            index = allocate_block(tag_bits, index_bits);
            if(lower_memory) lower_memory->read(address);
            memory[index] = tag_bits;
            // state_matrix[index] = VALID;
            if(replacement_policy == FIFO) {
                priority_matrix[index] = current_priority;   // update only when first accessed/allocated
                if(debug_mode) std::cout << cache_level << " update FIFO\n";
            }

        } else {
            if(debug_mode) std::cout << cache_level << " hit\n";
            // n_reads += 1;
        }

        if(state_matrix[index] != DIRTY) state_matrix[index] = VALID;
        // update priority_matrix
        if(replacement_policy == LRU) {
            priority_matrix[index] = current_priority;
            if(debug_mode) std::cout << cache_level << " update LRU\n";
        }
        if (debug_mode && replacement_policy == OPTIMAL) 
            std::cout << cache_level << " update optimal\n";
        current_priority = current_priority + 1;      // (current_priority+1) % MAX_PRIORITY ;

        return index; 
    }
        
    // Function to write_back a block
    void write_back(uint32_t address) {
        // if(!size || !set) return -1; // error code = -1 : cache doesnot exist
        
        if(debug_mode) std::cout << cache_level << " write : " << uint32ToHexChar(address) ; 
        n_writes += 1;

        uint32_t tag_bits = get_tag(address, indexBits, offsetBits);
        uint32_t index_bits = get_index(address, indexBits, offsetBits);
        
        if(debug_mode) std::cout << " (tag " << uint32ToHexChar(tag_bits) << ", index " << index_bits << ")\n" ; 

        // search the address 
        int index = search_memory(tag_bits, index_bits); 
        
        if(index == -1) {   // MISS. fetch from lower
            if(debug_mode) std::cout << cache_level << " miss\n";
            n_write_misses += 1;
            index = allocate_block(tag_bits, index_bits);
            if(lower_memory) {
                lower_memory->read(address);
            }
            
            memory[index] = tag_bits;
            if(replacement_policy == FIFO) {
                priority_matrix[index] = current_priority;   // update only when first accessed/allocated
                if(debug_mode) std::cout << cache_level << " update FIFO\n";
            }

        } else  {
            if(debug_mode) std::cout << cache_level << " hit\n";
            // n_writes += 1;
        }
        
        // update priority_matrix
        if(replacement_policy == LRU) {
            priority_matrix[index] = current_priority;
            if(debug_mode) std::cout << cache_level << " update LRU\n";
        }
        if (debug_mode && replacement_policy == OPTIMAL) 
            std::cout << cache_level << " update optimal\n";

        state_matrix[index] = DIRTY;    // written to this cache. but not consistant with the lower
        if(debug_mode) std::cout << cache_level << " set dirty\n";
        current_priority = current_priority + 1;         // (current_priority+1) % MAX_PRIORITY ;
    }

    double calculate_miss_rate () {
        if(cache_level == "L2" && n_reads) 
            return (double)n_read_misses/n_reads ;

        else if (cache_level == "L1" && (n_reads + n_writes)) 
            return (double)(n_read_misses + n_write_misses)/(n_reads + n_writes);

        return 0;
    }

    int get_traffic() {
        return n_read_misses+n_write_misses+n_write_backs ;
    }

    void print_stats(char bullet) {
        // if(size == 0) return;
        std::cout << (char)bullet << ". number of " << cache_level << " reads:\t" << n_reads << "\n";
        std::cout << (char)(bullet+1) << ". number of " << cache_level << " read misses:\t" << n_read_misses << "\n";
        std::cout << (char)(bullet+2) << ". number of " << cache_level << " writes:\t" << n_writes << "\n";
        std::cout << (char)(bullet+3) << ". number of " << cache_level << " write misses:\t" << n_write_misses << "\n";
        std::cout << (char)(bullet+4) << ". " << cache_level << " miss rate:\t" << calculate_miss_rate() << "\n";
        std::cout << (char)(bullet+5) << ". number of " << cache_level << " writebacks:\t" << n_write_backs << "\n";

    }

    // Display cache configuration
    void display_configuration() {
        std::cout << cache_level << " Cache Configuration:" << std::endl;
        std::cout << "  Cache Size: " << size << " bytes" << std::endl;
        std::cout << "  Associativity: " << assoc << std::endl;
        std::cout << "  Block Size: " << block_size << " bytes" << std::endl;
        std::cout << "  Replacement Policy: " << get_policy_name(replacement_policy) << std::endl;
        std::cout << "  Inclusion Property: " << get_inclusion_name(inclusion_property) << std::endl;
        // std::cout << "Done printing" << std::endl;
    }

    
    // Display cache contents
    void display_contents() {
        if(!size) return;
        int indx = 0;
        std::cout << "===== " << cache_level << " contents =====" << std::endl;
        for(int i=0; i<set; i++) {
            std::cout << "set \t" << i << ":\t" ;
            for (int j = 0; j<assoc; j++) {
                indx = i*assoc + j;
                std::cout << uint32ToHexChar(memory[indx]);
                if(state_matrix[indx] == DIRTY) std::cout << " D\t";
                else std::cout << "\t\t";
            }
            std::cout << std::endl;
        }


        // bool print_stats = false;

        // if(print_stats) {
        //     std::cout << "reads : " << n_reads << ", writes : " << n_writes << ", write_backs : " << n_write_backs;
        //     std::cout << "\nread_misses : " << n_read_misses << ", write_misses : " << n_write_misses << "\n";
        // }

    }

};


void display_simulation_conf(int block_size, int L1_size, int L1_assoc, int L2_size, int L2_assoc, 
    ReplacementPolicy replacement_policy, InclusionProperty inclusion_property, std::string trace_file) {
        
    std::cout << "===== Simulator configuration =====" << std::endl;
    std::cout << "BLOCKSIZE:\t\t" << block_size << std::endl;
    std::cout << "L1_SIZE:\t\t" << L1_size << std::endl;
    std::cout << "L1_ASSOC:\t\t" << L1_assoc << std::endl;
    std::cout << "L2_SIZE:\t\t" << L2_size << std::endl;
    std::cout << "L2_ASSOC:\t\t" << L2_assoc << std::endl;
    std::cout << "REPLACEMENT POLICY:\t" << get_policy_name(replacement_policy) << std::endl;
    std::cout << "INCLUSION PROPERTY:\t" << get_inclusion_name(inclusion_property) << std::endl;
    std::cout << "trace_file:\t\t" << trace_file << std::endl;
}


double calculate_AAT_L1(CACHE L1, double HT, int penalty) {
    if(L1.get_size() == 0) return 0;

    return HT + penalty * L1.calculate_miss_rate() ;
}


double calculate_AAT_L2(CACHE L1, CACHE L2, double HT_1, double HT_2, int penalty) {    // penalty is same for both in our case. So using one var
    if(L1.get_size() == 0) return 0;
    if(L2.get_size() == 0) return calculate_AAT_L1(L1, HT_1, penalty);

    return HT_1 + L1.calculate_miss_rate() * calculate_AAT_L1(L2, HT_2, penalty);
}


int main(int argc, char* argv[]) {
    if (graph && argc < 11) {
        std::cout << "Graph Usage: sim_cache <BLOCKSIZE> <L1_SIZE> <L1_ASSOC> <L2_SIZE> <L2_ASSOC> <REPLACEMENT_POLICY> <INCLUSION_PROPERTY> <trace_file> <L1_Hit_Time> <L2_Hit_Time>" << std::endl;
        return 1;

    } else if (argc < 9) {
        std::cout << "Simulation Usage: sim_cache <BLOCKSIZE> <L1_SIZE> <L1_ASSOC> <L2_SIZE> <L2_ASSOC> <REPLACEMENT_POLICY> <INCLUSION_PROPERTY> <trace_file>" << std::endl;
        return 1;
    }

    std::ofstream outFile;
    std::streambuf* originalCoutBuffer;
    if(output_to_file){
        outFile.open("output.txt");

        // Redirect std::cout to outFile
        originalCoutBuffer = std::cout.rdbuf(); // Save the original buffer
        std::cout.rdbuf(outFile.rdbuf()); // Redirect std::cout to outFile
    }



    int block_size = std::stoi(argv[1]);
    int L1_size = std::stoi(argv[2]);
    int L1_assoc = std::stoi(argv[3]);
    int L2_size = std::stoi(argv[4]) ;  // ? std::stoi(argv[4]) : 0 ;
    int L2_assoc = std::stoi(argv[5]) ; //  ? std::stoi(argv[5]) : 0 ;
    ReplacementPolicy policy = static_cast<ReplacementPolicy>(std::stoi(argv[6]));
    InclusionProperty inclusion = static_cast<InclusionProperty>(std::stoi(argv[7]));
    std::string trace_file = argv[8];

    
    double Hit_time_L1 = 0;       // std::stod(argv[9]);    // ns
    double Hit_time_L2 = 0;       // std::stod(argv[10]);   // ns
    int miss_penalty = 100; // ns
    if(graph) { 
        Hit_time_L1 = std::stod(argv[9]);
        Hit_time_L2 = std::stod(argv[10]);
    }

    if(simulation) 
        display_simulation_conf(block_size, L1_size, L1_assoc, L2_size, L2_assoc, policy, inclusion, trace_file);
 
    // Create L1 cache
    CACHE L1_cache("L1", L1_size, L1_assoc, block_size, policy, inclusion);

    CACHE L2_cache("L2", L2_size, L2_assoc, block_size, policy, inclusion);

    if(L2_size != 0) {  // if  L2 size is zero, it has no lower cache, thus no hierarchy
        L1_cache.set_memory_hierarchy(NULL, &L2_cache);
        // L1_cache.display_configuration();

        L2_cache.set_memory_hierarchy(&L1_cache, NULL);
        // L2_cache.display_configuration();
    }

    
    int set = L1_size / (L1_assoc * block_size);
    int offsetBits = std::log2(block_size);
    int indexBits = std::log2(set);


    // for Optimal Replacement Policy
    if (policy == OPTIMAL) {

        std::ifstream file(trace_file);
        if (!file.is_open()) {
            std::cout << "Failed to open trace file : " << trace_file << "\tfor optimal \n";
            return 1;
        }

        int j = 0;
        std::string line_o;
        while (std::getline(file, line_o)) {
            std::istringstream iss_o(line_o);
            char op;            // Operation: 'r' for read, 'w' for write_back
            std::string address_hex;

            if (!(iss_o >> op >> address_hex)) {
                continue;       // Skip invalid lines
            }

            uint32_t address = hexAddrToUint32(address_hex);      
            uint32_t addr_key = get_address(get_tag(address, indexBits, offsetBits),
                                            get_index(address, indexBits, offsetBits), 
                                            0, indexBits, offsetBits); 

            // std::cout << address_hex << "\t:\t" << uint32ToHexChar_wzero(addr_key) << "\n";
            // address = address << offsetBits ;
            // address = address >> offsetBits ;

            access_cycle_tracker[addr_key].push_back(j);
            j++;
            // if(j == 100) break;      // exit(0);
        }
        
        file.close();

        // int infinite_val = j+99999;    // will use this as infinite

        for (auto it_map = access_cycle_tracker.begin(); it_map != access_cycle_tracker.end(); ++it_map) {
            it_map->second.push_back(INFINITE);
            // if (get_index(it_map->first, indexBits, offsetBits) == 31) {
            //     std::cout << "index 31 : " << uint32ToHexChar(get_tag(it_map->first, indexBits, offsetBits)) << "\t:\t" ;
            //     for (auto it_vect = it_map->second.begin(); it_vect != it_map->second.end(); ++it_vect) {
            //         std::cout << *it_vect << ", " ;
            //     }
            //     std::cout << "\n\n";
            // }
        }
    }



    // start reading trace file for input and main process
    std::ifstream file(trace_file);
    if (!file.is_open()) {
        std::cout << "Failed to open trace file: " << trace_file << std::endl;
        return 1;
    }

    int i = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char op;            // Operation: 'r' for read, 'w' for write_back
        std::string address_hex;

        if (!(iss >> op >> address_hex)) {
            continue;       // Skip invalid lines
        }

        uint32_t whole_address = hexAddrToUint32(address_hex);        
        uint32_t address = get_address(get_tag(whole_address, indexBits, offsetBits),
                                            get_index(whole_address, indexBits, offsetBits), 
                                            0, indexBits, offsetBits);

        if(debug_mode) std::cout << "----------------------------------------\n" ;

        if (op == 'r') {
            if(debug_mode) std::cout << "# " << i+1 <<" : read " << address_hex << "\n";
            L1_cache.read(address) ;  
        } else if (op == 'w') {
            if(debug_mode) std::cout << "# " << i+1 <<" : write " << address_hex << "\n";
            L1_cache.write_back(address) ;
        }
        
        if(debug_mode && i == max_debug_line) break;
        i++;
    }

    file.close();
    
    if(simulation) {
        L1_cache.display_contents();
        L2_cache.display_contents();
    }

    if(graph) {
        std::cout << L1_cache.calculate_miss_rate() << "\n";
        if(L2_size == 0) 
            std::cout << calculate_AAT_L1(L1_cache, Hit_time_L1, miss_penalty) << "\n";
        else 
            std::cout << calculate_AAT_L2(L1_cache, L2_cache, Hit_time_L1, Hit_time_L2, miss_penalty) << "\n";
    }
    
    if(simulation) {
        // print simulations
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "===== Simulation results (raw) =====\n";
        L1_cache.print_stats('a');
        L2_cache.print_stats('g');

        std::cout << "m. total memory traffic:\t" ;
        if(L2_size) std::cout << L2_cache.get_traffic() + L1_cache.get_write_backs_to_main() << "\n";
        else std::cout << L1_cache.get_traffic() << "\n"; 
    }

    L1_cache.free_memories();
    L2_cache.free_memories();


    if(output_to_file){
        // Restore the original std::cout buffer (console output)
        std::cout.rdbuf(originalCoutBuffer);        
        // Close the file
        outFile.close();
    }
    
    // std::cout << "Done Processing " << i << " requests\n";

    return 0;
}
