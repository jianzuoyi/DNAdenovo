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
#include "gzstream.h"

using namespace std;

//#define KMERLEN 34;
//#define KMERFREQ 30;
#define KMERLEN 56;
#define KMERFREQ 8;

//define the entity data structure (elements of the hashset-array)
typedef struct Entity{
	uint64_t kmer:KMERLEN;
	uint64_t freq:KMERFREQ;
};


//define the hashSet data structure(include data and control parameters)
typedef struct {
	uint32_t e_size; //the byte size for each entity 
	uint64_t size; //the space size of the array
	uint64_t count; //the number of current stored entities
	uint64_t count_conflict; //the number of current stored entities
	uint64_t max; //the max number of entities allowed to store
	float load_factor; //the ratio (allowed entity number / total array size) 
	uint64_t iter_ptr; //the pointer to the current processing entity 
	Entity *array;  //the arrary to store the body data
	uint8_t *nul_flag; //the array to store null flags, one bit for an entity
	uint8_t *del_flag; //the array to store delete flags, one bit for an entity
} HashSet;


uint64_t HashCode(uint64_t key);
uint64_t hash_code(Entity*& entity);
int hash_equal(Entity*& a, Entity* b);
int free_hash(HashSet *set);

int is_prime(uint64_t num);
uint64_t find_next_prime(uint64_t num);
HashSet *init_hashset(uint64_t init_size, float load_factor);

int is_entity_null(uint8_t *nul_flag, uint64_t idx);
int set_entity_fill(uint8_t *nul_flag, uint64_t idx);
int is_entity_delete(uint8_t *del_flag, uint64_t idx);
int set_entity_delete(uint8_t *del_flag, uint64_t idx);

void enlarge_hashset (HashSet*& set, int num);
int add_hashset(HashSet*& set, Entity*& entity);
uint64_t get_hashset(HashSet*& set, Entity*& entity);
int exists_hashset (HashSet*& set, Entity*& entity);
int delete_hashset (HashSet*& set, Entity*& entity);
//void print_hashset(HashSet *set);
uint64_t print_hashset(HashSet *set, ogzstream& fp, uint64_t*& freqy);

uint64_t get_hashset(HashSet*& set, Entity*& entity);


#endif


