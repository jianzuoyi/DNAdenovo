#include "hashSet.h"

//Jenkins' hash function for 64-bit integers, convert Kmer into 64-bit integers
uint64_t hash_code(Entity*& entity){
	uint64_t key=entity->kmer;
        key += ~(key << 32);
        key ^= (key >> 22);
        key += ~(key << 13);
        key ^= (key >> 8);
        key += (key << 3);
        key ^= (key >> 15);
        key += ~(key << 27);
        key ^= (key >> 31);
	return key;
}

uint64_t HashCode(uint64_t key){
        key += ~(key << 32);
        key ^= (key >> 22);
        key += ~(key << 13);
        key ^= (key >> 8);
        key += (key << 3);
        key ^= (key >> 15);
        key += ~(key << 27);
        key ^= (key >> 31);
	return key;
}


//The hash key comparison function 
int hash_equal(Entity*& a, Entity* b)
{ 
	return (a->kmer == b->kmer); 
}

//free the memory of HashSet data structure
int free_hash(HashSet *set)
{ 
	free(set->array);
	free(set->nul_flag);
	free(set->del_flag);
	free(set); 
	return 1; 
}

//judge whether a number is a prime
int is_prime(uint64_t num)
{
	uint64_t i, max;
	if (num < 4) return 1;
	if (num % 2 == 0) return 0;
	max = (uint64_t)sqrt((float)num);
	for (i=3;i<max;i+=2){ if (num % i == 0) return 0; }
	return 1;
}

//find the next prime number for a given number
uint64_t find_next_prime(uint64_t num)
{
	if (num % 2 == 0) num ++;
	while (1) { if (is_prime(num)) return num; num += 2; }
}


//Initialization of the HashSet data structure
HashSet *init_hashset(uint64_t init_size, float load_factor)
{
	HashSet *set = (HashSet*) malloc(sizeof(HashSet));
	
	if (init_size < 3) init_size = 3;
	else init_size = find_next_prime(init_size);

	set->e_size = sizeof(Entity);
	set->size   = init_size;
	set->count  = 0;
	set->count_conflict = 0;

	if (load_factor <= 0) load_factor = 0.25f;
	else if (load_factor >= 1) load_factor = 0.75f;	
	
	set->load_factor = load_factor;
	set->max    	 = (uint64_t) (set->size * load_factor);
	set->iter_ptr    = 0;

	set->array = (Entity*) malloc(set->size * set->e_size); 
	if (set->array == NULL) { cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << "\n"; exit(1); }
	memset(set->array, 0, set->size*set->e_size);

	set->nul_flag = (uint8_t*) malloc(set->size / 8 + 1); 
	if (set->nul_flag == NULL) { cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << "\n"; exit(1); }
	memset(set->nul_flag,0,set->size / 8 + 1);  //initialization

	set->del_flag = (uint8_t*) malloc(set->size / 8 + 1); 
	if (set->del_flag == NULL) { cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << "\n"; exit(1); }
	memset(set->del_flag,0,set->size / 8 + 1);  //initialization

	return set;
}

//check whether the positon of 'idx' in array is empty
int is_entity_null(uint8_t *nul_flag, uint64_t idx) 
{
	int value = ( nul_flag[idx/8] >> (7-idx%8) ) & 0x1u;
	return (1-value);

}

//set the nul_flag to 1 when an entity is put in the position of array
int set_entity_fill(uint8_t *nul_flag, uint64_t idx) 
{
	uint8_t a[8] = {128,64,32,16,8,4,2,1};
	nul_flag[idx/8] |= a[idx%8];
}

//check whether the entity in the position of 'idx' has been deleted
int is_entity_delete(uint8_t *del_flag, uint64_t idx) 
{
	int value = (del_flag[idx/8] >> (7-idx%8)) & 0x1u;
	return value;
}

//set the del_flag to 1 when an entity was deleted from position of idx in array
int set_entity_delete(uint8_t *del_flag, uint64_t idx) 
{
	uint8_t a[8] = {128,64,32,16,8,4,2,1};
	del_flag[idx/8] |= a[idx%8];
}

//enlarge the HashSet memory space (*2 each time) 
void enlarge_hashset (HashSet*& set, int num) 
{
	uint64_t old_size = set->size;
	uint64_t new_size = set->size;
	do{ new_size = find_next_prime(new_size * 2); } while(new_size * set->load_factor < set->count + num);
	
	set->size  = new_size;
	set->max   = (uint64_t)(new_size * set->load_factor+0.5);
	set->array = (Entity*) realloc(set->array, new_size*set->e_size);
	if (set->array == NULL)
	{
		cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << ".\n";
		exit(1);
	}

	uint8_t *nul_flag, *del_flag;
	nul_flag = set->nul_flag;
	del_flag = set->del_flag;

	set->nul_flag = (uint8_t*) malloc(new_size/8 + 1);
	if (set->nul_flag == NULL)
	{
		cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << ".\n";
		exit(1); 
	}
	memset(set->nul_flag, 0, new_size/8+1); //initialization
	set->del_flag = (uint8_t*) malloc(new_size/8 + 1);
	if (set->del_flag == NULL)
	{
		cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << ".\n";
		exit(1); 
	}
	memset(set->del_flag, 0, new_size/8+1); //initialization

	Entity *tmp1, *tmp2;
	tmp1 = (Entity*) malloc(set->e_size);
	tmp2 = (Entity*) malloc(set->e_size);
	if (tmp1 == NULL || tmp2 == NULL)
	{
		cerr << "Out of memory! Line: " << __LINE__ << ". Function: " << __FUNCTION__ << ".\n";
		exit(1); 
	}

	for (uint64_t i=0; i<old_size; i++)
	{
		if ((is_entity_null(nul_flag, i)==1) || (is_entity_delete(del_flag, i)==1)) continue;
		memcpy(tmp1, set->array+i, set->e_size);
		set_entity_delete(del_flag, i);
		while(1)
		{
			uint64_t hc = hash_code(tmp1) % set->size;
			while (is_entity_null(set->nul_flag, hc) == 0) { hc = (hc + 1) % set->size; }
			set_entity_fill(set->nul_flag, hc);
			if ((hc < old_size) && (is_entity_null(nul_flag, hc)==0) && (is_entity_delete(del_flag, hc)==0))
			{/* the position was occupied by an old entity, replace it with new entity, 
					then find a new positon for the old entity */
				memcpy(tmp2, set->array+hc, set->e_size);
				memcpy(set->array+hc, tmp1, set->e_size);
				memcpy(tmp1, tmp2, set->e_size);
				set_entity_delete(del_flag, hc);
			}
			else
			{/* found an empty position */
				memcpy(set->array+hc, tmp1, set->e_size);
				break;
			}
		}
	}

	free(nul_flag);
	free(del_flag);
	free(tmp1);
	free(tmp2);
}


//Find a position(hc,hc+1,hc+2,...), and put the entity there.
//When the key is not existed, insert it in an empty(nul_flag is 0) position.
//When the key is already existed, just update the value of entity, but do not change the key  
int add_hashset(HashSet*& set, Entity*& entity)
{
	if (set->count+1 > set->max)
	{
		enlarge_hashset(set, 1);
	}

	uint64_t hc = hash_code(entity) % set->size;

	do{
		if (is_entity_null(set->nul_flag, hc))
		{
			memcpy(set->array+hc, entity, set->e_size);
			set_entity_fill(set->nul_flag, hc);
			set->count++;
			return 1;
		}
		else
		{
			if (hash_equal(entity,set->array+hc))
			{
				set->array[hc].freq++;
				return 1;
			}
		}
		
		set->count_conflict ++;
		if (hc + 1 == set->size) 
		{ 
			hc = 0; 
		}
		else 
		{
			hc = hc + 1;
		}
	} while(1);

	
	return 0;
}


//return the array index, if the entity is existed 
//return the array size, if the entity is not existed 
uint64_t get_hashset(HashSet*& set, Entity*& entity)
{
	uint64_t hc = hash_code(entity) % set->size;
	do{
		if (is_entity_null(set->nul_flag, hc) == 1)
		{
			return set->size;
		}
		else if (hash_equal(entity,set->array+hc))
		{
			if (is_entity_delete(set->del_flag, hc) == 0)
			{
				return hc;
			}
			else
			{
				return set->size;
			}
		}

		if (hc + 1 == set->size) { hc = 0; }
		else hc = hc + 1;
	} while(1);
	return 0;
}


//check whether an entity existes in the HashSet
int exists_hashset (HashSet *set, Entity *entity)
{
	uint64_t idx = get_hashset(set, entity);
	int is_exists = 0;
	if (idx != set->size) is_exists = 1;
	return is_exists;
}

//delete an entity from the HashSet, by setting the del_flag to 1
int delete_hashset (HashSet * &set, Entity * &entity)
{
	uint64_t idx = get_hashset(set, entity);
	if (idx != set->size)
	{
		set_entity_delete(set->del_flag, idx);
		set->count--;
	}
}

//print kmer sequence and value
void print_hashset(HashSet *set)
{
	Entity *array = set->array;
	for (uint64_t i=0; i<set->size; i++)
	{	
		if (!is_entity_null(set->nul_flag, i))
		{	
			cout << set->array[i].kmer << "\t" << set->array[i].freq << "\n";
		}
	}
	
}
