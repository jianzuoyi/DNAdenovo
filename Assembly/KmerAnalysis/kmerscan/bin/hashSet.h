#ifndef __HASHSET_H_
#define __HASHSET_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
using namespace std;

//define the entity data structure (elements of the hashset-array)
typedef struct {
	uint64_t low;
	uint64_t high:52;	//store Kmer-number if Kmer>32.
	uint64_t freq:4;	//store frequency
	uint64_t leftBase:3;	//base on the left side of this kmer
	uint64_t rightBase:3;	//base on the right side of thisi kmer
	uint64_t leftBranch:1;	//if this kmer has more than one kind of base on its left side, leftBranch=1
	uint64_t rightBranch:1;	//if this kmer has more than one kind of base on its right side, rightBranch=1
}Entity;

//define the hashSet data structure(include data and control parameters)
typedef struct {
	uint32_t e_size; //the byte size for each entity 
	Entity *array;  //the array to store the body data
	uint64_t size; //the size of the array declared above
	uint64_t count; //the number of current stored entities
	uint64_t count_conflict; //the number of conflict
	uint64_t max; //the max number of entities allowed to store
	uint64_t iter_ptr; //the pointer to the current processing entity 
	float load_factor; //the ratio (allowed entity number / total array size) 
	uint8_t *nul_flag; //the array to store null flags, one bit for an entity
	uint8_t *del_flag; //the array to store delete flags, one bit for an entity
} HashSet;


uint64_t hash_code(Entity *e);
int hash_equal(Entity *a, Entity *b);
int free_hash(HashSet *set);

int is_prime(uint64_t num);
uint64_t find_next_prime(uint64_t num);
HashSet *init_hashset(uint64_t init_size, float load_factor);

int is_entity_null(uint8_t *nul_flag, uint64_t idx);
int set_entity_fill(uint8_t *nul_flag, uint64_t idx);
int is_entity_delete(uint8_t *del_flag, uint64_t idx);
int set_entity_delete(uint8_t *del_flag, uint64_t idx);

void enlarge_hashset (HashSet *set, uint64_t num);
int add_hashset(HashSet *set, Entity *entity);
uint64_t get_hashset(HashSet *set, Entity *entity);
int exists_hashset (HashSet *set, Entity *entity);
int delete_hashset (HashSet *set, Entity *entity);
void print_hashset(HashSet *set, ofstream& kmer_outfile, ofstream &kmerFreq_outfile);
uint64_t get_freq(HashSet* set, uint64_t* kmer, uint64_t *leftBranch, uint64_t *rightBranch);


#endif


