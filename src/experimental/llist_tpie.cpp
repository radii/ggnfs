/**************************************************************
 * llist_tpie.cpp
 * Routines for efficient operation on (large) arrays with a
 * variable number of s32s in each entry.
 * Copyright 2005 Anton Korobeynikov.
 ***************************************************************/

/*  This file is part of GGNFS.
 *
 *   GGNFS is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   GGNFS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with GGNFS; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "llist_tpie.h"
#include "templates.h"

template <class BTECOLL>
_LBListField<BTECOLL>::_LBListField(AMI_collection_single<BTECOLL>* pcoll, AMI_bid ID, entry_count_t maxEntries):
	AMI_block<s32, _LBList_entry_info, BTECOLL>(pcoll, 0, ID) {

	if (ID == 0) {
		info()->entries = 0;
	}

	_maxEntries = el.capacity();
	if (maxEntries <= _maxEntries && maxEntries != 0) {
		_maxEntries = maxEntries;
	}
}


template <class BTECOLL>
_LBList<BTECOLL>::_LBList(size_t lbf, persistence per):_entries(NULL), _count(0), _maxEntries(0) {
	char *base_name = tpie_tempnam("LBList");
	_collection=new AMI_collection_single<BTECOLL>(base_name,AMI_READ_WRITE_COLLECTION,lbf);
	_collection->persist(per);
	_entry_cache = new entry_cache_t(200000,8);
	printf("Block size: %d\n",_collection->block_size());
};

template <class BTECOLL>
_LBList<BTECOLL>::_LBList(llist_t* L, size_t lbf, persistence per):_entries(NULL), _count(0), _maxEntries(0) {
	char *base_name = tpie_tempnam("LBList");
	_collection=new AMI_collection_single<BTECOLL>(base_name,AMI_READ_WRITE_COLLECTION,lbf);
	_collection->persist(per);
	_entry_cache = new entry_cache_t(200000,8);
	printf("Block size: %d\n",_collection->block_size());

	for (s32 i=0; i<L->numFields; ++i) {
		reallocateEntries();
		_LBListField<BTECOLL>* newEntry = new _LBListField<BTECOLL> (_collection,(AMI_bid)0);
		_entries[_count] = newEntry->bid();
		for (int j=L->index[i]; i<L->index[i+1]; ++i) {
			(newEntry->el)[j-L->index[i]]=L->data[j];
		}
		newEntry->info()->entries=L->index[i+1]-L->index[i];
		newEntry->persist(PERSIST_PERSISTENT);
		++_count;
		delete newEntry;
	}	
}

template <class BTECOLL>
_LBList<BTECOLL>::~_LBList() {
	delete[] _entries;
	_entry_cache->flush();
	delete _entry_cache;
	delete _collection;
}

template <class BTECOLL>
void _LBList<BTECOLL>::clearList() {
	_entry_cache->flush();
	for (u32 i=0; i<_count; ++i) {
		_LBListField<BTECOLL>* cField= new _LBListField<BTECOLL>(_collection,_entries[i],0);
		cField->persist(PERSIST_DELETE);
		delete cField;
	}
	delete[] _entries;
}

template <class BTECOLL>
_LBListField<BTECOLL>* _LBList<BTECOLL>::fetchField (unsigned int index) {
	if (index >= _count) {
		std::cout << "Error in LList::fetchField: index out of range: " << index <<" out of " << _count << endl ;
		assert(false);
	}
	_LBListField<BTECOLL>* newEntry;

	if (_entries[index]==0) {
		cout << "Error in LList::fetchField: ID == 0. Index=" << index << endl;
		assert(false);
	}
	if (!_entry_cache->read(_entries[index], newEntry)) {
		newEntry = new _LBListField<BTECOLL>(_collection,_entries[index],0);
	}

	return newEntry;
}

template <class BTECOLL>
void _LBList<BTECOLL>::releaseField (_LBListField<BTECOLL>* entry) {
	if (entry->persist() == PERSIST_DELETE) {
		delete entry;
	} else {
		_entry_cache->write(entry->bid(), entry);
	}
}

template <class BTECOLL>
void _LBList<BTECOLL>::reallocateEntries() {
	if ((_count + 1) > _maxEntries) {
		printf("reallocateEntries()::_maxEntries=%d\n",_maxEntries);
		_maxEntries+=BLOCK_SIZE;
		AMI_bid* newentries = new AMI_bid[_maxEntries];
		memcpy(newentries,_entries,_count*sizeof(AMI_bid));
		delete[] _entries;
		_entries = newentries;
	}
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::appendField(s32 *entries, s32 numEntries) {
	AMI_err result = AMI_ERROR_NO_ERROR;

	reallocateEntries();
	_LBListField<BTECOLL>* newEntry = new _LBListField<BTECOLL> (_collection,(AMI_bid)0);
	_entries[_count] = newEntry->bid();
	for (int i=0; i<numEntries; ++i) {
		(newEntry->el)[i]=entries[i];
	}
	newEntry->info()->entries=numEntries;
	newEntry->persist(PERSIST_PERSISTENT);
	++_count;
	releaseField(newEntry);
	return result;
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::appendField(const _LBListField<BTECOLL>& oldEntry) {
	AMI_err result = AMI_ERROR_NO_ERROR;

	reallocateEntries();

	_LBListField<BTECOLL>* cField = new _LBListField<BTECOLL>(_collection,(AMI_bid)0);

	_entries[_count] = cField->bid();
	for (u32 i=0; i<(oldEntry.info())->entries; ++i) {
		(cField->el)[i]=oldEntry.el[i];
	}
	cField->info()->entries=(oldEntry.info())->entries;
	cField->persist(PERSIST_PERSISTENT);
	++_count;
	releaseField(cField);
	return result;
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::appendToField(s32 entry, s32 *entries, s32 numEntries) {
	AMI_err result = AMI_ERROR_NO_ERROR;

	reallocateEntries();
	_LBListField<BTECOLL>* newEntry = fetchField(entry);
	for (int i=0; i<numEntries; ++i) {
		(newEntry->el).insert(entries[i],(newEntry->info()->entries)++);
	}
	releaseField(newEntry);
	return result;
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::deleteFields(s32 *fields, s32 numFields) {
	AMI_err result = AMI_ERROR_NO_ERROR;

	for (int i=0; i<numFields; ++i) {
		if ((u32)fields[i] >= _count) {
			std::cout << "Error in LList::deleteFields: index out of range: " << fields[i] <<" out of " << _count << endl ;
			assert(false);
		}
		_LBListField<BTECOLL>* deleteEntry = fetchField(fields[i]);
		deleteEntry->persist(PERSIST_DELETE);
		releaseField(deleteEntry);
	}

	s32 cIdx=fields[0];
	for (int i=0; i<numFields; ++i) {
		s32 moveEnd = (i<(numFields-1)) ? (fields[i+1]):(_count);
		for (int j=fields[i]+1; j<moveEnd; ++j, ++cIdx) {
			_entries[cIdx]=_entries[j];
		}
	}
	_count-=numFields;
	return result;
}

template <class BTECOLL>
int _LBList<BTECOLL>::numCommonEntries(s32 field0, s32 field1) {
	u32 i, j, num=0;

	assert (field0 < (s32)_count && "index out of range");
	assert (field1 < (s32)_count && "index out of range");

	_LBListField<BTECOLL>* entry0 = fetchField(field0);
	_LBListField<BTECOLL>* entry1 = fetchField(field1);

	for (i=0; i<entry0->info()->entries; i++) {
		for (j=0; j<entry1->info()->entries; j++) {
			if ((entry0->el)[i]==(entry1->el)[j]) {
				num++;
			}
		}
	}
	releaseField(entry0);
	releaseField(entry1);
	return num;
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::catFields(lpair_t *pairs, s32 numPairs, int mod2) {
	s32 newEntries[MAX_FIELD_ENTRIES], numNewEntries;
	AMI_err result = AMI_ERROR_NO_ERROR;

	for (s32 i=0; i<numPairs; ++i) {
		printTmp("LList::catFields: total %" PRId32 " pairs done out of %" PRId32 ".",i,numPairs);
		assert(pairs[i].x < (s32)_count && "index out of range!");
		assert(pairs[i].y < (s32)_count && "index out of range!");

		_LBListField<BTECOLL>* entry0 = fetchField(pairs[i].x);
		_LBListField<BTECOLL>* entry1 = fetchField(pairs[i].y);

		numNewEntries=0;
		for (u32 j=0; (j<entry0->info()->entries) && (numNewEntries < MAX_FIELD_ENTRIES); j++) {
			newEntries[numNewEntries++] = entry0->el[j];
		}
		for (u32 j=0; (j<entry1->info()->entries) && (numNewEntries < MAX_FIELD_ENTRIES); j++) {
			newEntries[numNewEntries++] = entry1->el[j];
		}
		if (numNewEntries >= MAX_FIELD_ENTRIES) {
			printf("MAX_FIELD_ENTRIES exceeded (i=%" PRId32 ", c0=%" PRId32 ", c1=%" PRId32 ".). Ignoring...\n",i,pairs[i].x,pairs[i].y);
			return AMI_ERROR_GENERIC_ERROR;
		}

		qsort(newEntries, numNewEntries, sizeof(s32), cmpS32s);
		int newEndLoc=0;
		int i=0;
		int j=0;
		while (i<numNewEntries) {
			/* How many copies of newEntries[i] are there? */
			j=i+1;
			while ((j<numNewEntries) && (newEntries[j]==newEntries[i])) {
				j++;
			}
			if (mod2) {
				if ((j-i)%2 == 0) {
					/* An even number of them, so omit them altogether. */
					i=j;
				} else {
					/* Copy one of them over. */
					newEntries[newEndLoc++] = newEntries[i];
					i=j;
				}
			} else {
				/* Copy one of them over. */
				newEntries[newEndLoc++] = newEntries[i];
				i=j;
			}
		}

		for (int i=0; i<newEndLoc; ++i) {
			(entry0->el)[i]=newEntries[i];
		}
		(entry0->info())->entries=newEndLoc;

		releaseField(entry0);
		releaseField(entry1);
	}
	return result;
}

template <class BTECOLL>
AMI_err _LBList<BTECOLL>::getSortOnFieldSize(s32 *fields) {
	s32 counts[MAX_LIST_ENTRIES];
	s32 c_idx;
	AMI_err result = AMI_ERROR_NO_ERROR;

	memset(counts,0,sizeof(s32)*MAX_LIST_ENTRIES);
	for (u32 i=0; i<_count; ++i) {
		c_idx=fieldSize(i);
		if (c_idx<MAX_LIST_ENTRIES) {
			++counts[c_idx+1];
		} else {
			printf("getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
			return AMI_ERROR_GENERIC_ERROR;
		}
	}
	for (u32 i=1; i<MAX_LIST_ENTRIES; ++i) {
		counts[i]+=counts[i-1];
	}

	for (u32 i=0; i<_count; ++i) {
		c_idx=fieldSize(i);
		if (c_idx<MAX_LIST_ENTRIES) {
			fields[counts[c_idx]++]=i;
		} else {
			printf("getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
			return AMI_ERROR_GENERIC_ERROR;
		}
	}
	return result;
}

// FIXME: rewrite ;)
template <class BTECOLL>
s32 _LBList<BTECOLL>::getMaxEntry() {
	s32 N=0;
	for (u32 i=0; i<_count; ++i) {
		_LBListField<BTECOLL>* cField = fetchField(i);
		for (u32 j=0; j<cField->info()->entries;++j) {
			if ((cField->el[j] > N) && (cField->el[j] != (s32)BAD_LP_INDEX)) {
				N=cField->el[j];
			}
		}
		releaseField(cField);
	}
	return N;
}

// FIXME: rewrite ;)
template <class BTECOLL>
s32 _LBList<BTECOLL>::Weight() {
	s32 N=0;
	for (u32 i=0; i<_count; ++i) {
		_LBListField<BTECOLL>* cField = fetchField(i);
		N+=cField->info()->entries;
		releaseField(cField);
	}
	return N;
}

// FIXME: rewrite ;)
template <class BTECOLL>
s32 _LBList<BTECOLL>::fieldSize(u32 field) {
	_LBListField<BTECOLL>* cField = fetchField(field);	
	s32 res = cField->info()->entries;
	releaseField(cField);
	return res;
}

template <class BTECOLL>
llist_t * _LBList<BTECOLL>::toLList_t() {
	llist_t* P = (llist_t *)lxmalloc(sizeof(llist_t),1);
	ll_init(P,Weight(),numFields());
	for (u32 i=0; i<_count; ++i) {
		// do bulk load. Without any cache.
		_LBListField<BTECOLL>* cField = fetchField(i);
		for (u32 j=0; j<cField->info()->entries;++j) {
			P->data[P->index[i]+j]=cField->el[j];
		}
		P->index[i+1]=P->index[i]+cField->info()->entries;
		releaseField(cField);
	}
	return P;	
}

template <class BTECOLL>
void _LBList<BTECOLL>::sortEntries(u32 field) {
	_LBListField<BTECOLL>* cField = fetchField(field);
	s32* entries = new s32[cField->info()->entries];
	for (u32 i=0; i<cField->info()->entries; ++i) {
		entries[i] = cField->el[i];
	}
	qsort(entries, cField->info()->entries, sizeof(s32), cmpS32s);
	for (u32 i=0; i<cField->info()->entries; ++i) {
		cField->el[i] = entries[i];
	}
	delete[] entries;
	releaseField(cField);
}

_LSListField::_LSListField(s32* entries, s32 numEntries):_m_entries(NULL) {
	_m_entries = new std::vector<s32>;
	for (int i=0; i<numEntries; ++i) {
		_m_entries->push_back(entries[i]);
	}
}

_LSListField::~_LSListField() {
	delete _m_entries;
}	

void _LSListField::appendEntry(s32 entry) {
	if (!_m_entries) {
		_m_entries = new std::vector<s32>;	
	}
	_m_entries->push_back(entry);
}

_LSList::_LSList(persistence per):_m_fields(NULL), _m_per(per) {
	char *base_name = tpie_tempnam("LSList");

	_m_fields = new std::vector<off_t>;
	_m_fields->push_back(0);
	_m_baseStream = new AMI_stream<s32>(base_name);
	_m_baseStream->persist(per);
}

_LSList::_LSList(llist_t* L, persistence per):_m_fields(NULL), _m_per(per) {
	char *base_name = tpie_tempnam("LSList");

	_m_fields = new std::vector<off_t>;
	_m_fields->push_back(0);
	_m_baseStream = new AMI_stream<s32>(base_name);
	_m_baseStream->persist(per);

	for (s32 i=0; i<L->numFields; ++i) {
		appendField(&(L->data[L->index[i]]),L->index[i+1]-L->index[i]);
	}	

}

_LSList::~_LSList() {
	delete _m_fields;
	delete _m_baseStream;
}


_LSListField* _LSList::fetchField(u32 index) {
	// FIXME: assert
	TPIE_OS_OFFSET fs = fieldSize(index);
	s32* entries = new s32[fs];
	memset(entries,0,fs*sizeof(s32));
	_m_baseStream->seek(_m_fields->at(index));
	_m_baseStream->read_array(entries,&fs);
	_LSListField* F = new _LSListField(entries,fs);
	delete[] entries;
	return F;
}

void _LSList::releaseField(_LSListField* field) {
	// FIXME: assert
	delete field;
}

AMI_err	_LSList::appendField(s32 *entries, s32 numEntries) {
	AMI_err result = AMI_ERROR_NO_ERROR;

	u32 streamLen = _m_baseStream->stream_len();
	_m_baseStream->seek(streamLen);
	_m_baseStream->write_array(entries,numEntries);
	u32 cSize = Weight();
	_m_fields->push_back(cSize+numEntries);

	return result;
}

AMI_err	_LSList::deleteFields(s32 *fields, s32 numFields) {
	u32 buffSize = 0;
	s32* buff = NULL;

	char *base_name = tpie_tempnam("LSList");
	AMI_err result = AMI_ERROR_NO_ERROR;

	AMI_stream<s32>* newBaseStream = new AMI_stream<s32>(base_name);
	newBaseStream->persist(_m_per);

	s32 cdfield = 0;
	_m_baseStream->seek(0);
	for (u32 i=0; i<_m_count();++i) {
		TPIE_OS_OFFSET fs = fieldSize(i);
		if (fs>buffSize) {
			buffSize = fs;
			delete[] buff;
			buff = new s32[buffSize];
			memset(buff,0,buffSize);
		}
		_m_baseStream->read_array(buff,&fs);
		if (cdfield>=numFields || (s32)i != fields[cdfield]) {
			newBaseStream->write_array(buff,fs);
			_m_fields->at(i-cdfield+1)=_m_fields->at(i-cdfield)+fs;
		} else {
			++cdfield;
		}
	}
	_m_fields->erase(_m_fields->begin()+(_m_count()-numFields+1),_m_fields->end());
	delete[] buff;	delete _m_baseStream;
	_m_baseStream = newBaseStream;

	return result;
}


int _LSList::numCommonEntries(s32 field0, s32 field1) {
	u32 i, j, num=0;

	_LSListField* entry0 = fetchField(field0);
	_LSListField* entry1 = fetchField(field1);

	for (i=0; i<entry0->numEntries(); i++) {
		for (j=0; j<entry1->numEntries(); j++) {
			if (entry0->at(i)==entry1->at(j)) {
				++num;
			}
		}
	}
	releaseField(entry0);
	releaseField(entry1);
	return num;
}

inline bool checkPair(s32* pairs_x, s32 count, s32 y) {
	return (bsearch(&y,pairs_x,count,sizeof(s32),cmpS32s)==NULL);
}

AMI_err	_LSList::catFields(lpair_t *pairs, s32 numPairs, bool mod2) {
	s32 newEntries[MAX_FIELD_ENTRIES], numNewEntries;
	AMI_err result = AMI_ERROR_NO_ERROR;
	char *base_name = tpie_tempnam("LSList");
	u32 buffSize = 0;
	s32* buff = NULL;
	
	std::vector<off_t>* newOffsets = new std::vector<off_t>;
	newOffsets->push_back(0);

	
	AMI_stream<s32>* newBaseStream = new AMI_stream<s32>(base_name);
	newBaseStream->persist(_m_per);

	_LSListField** pairs_y = new _LSListField*[numPairs];
	s32* pairs_x = new s32[numPairs];
	
	for (s32 i=0; i<numPairs; ++i) {
		pairs_y[i] = fetchField(pairs[i].y);
		pairs_x[i] = pairs[i].x;
	}

	_m_baseStream->seek(0);
	s32 cdfield = 0;
	for (u32 i=0; i<_m_count(); ++i) {
		TPIE_OS_OFFSET fs = fieldSize(i);
		if (fs>buffSize) {
			buffSize = fs;
			delete[] buff;
			buff = new s32[buffSize];
		}
		_m_baseStream->read_array(buff,&fs);
		if (cdfield<numPairs && (s32)i==pairs[cdfield].x && checkPair(pairs_x,cdfield,pairs[cdfield].y)) {
			numNewEntries=0;
			for (u32 j=0; (j<fs) && (numNewEntries < MAX_FIELD_ENTRIES); j++) {
				newEntries[numNewEntries++] = buff[j];
			}
			for (u32 j=0; (j<pairs_y[cdfield]->numEntries()) && (numNewEntries < MAX_FIELD_ENTRIES); j++) {
				newEntries[numNewEntries++] = pairs_y[cdfield]->at(j);
			}
			if (numNewEntries >= MAX_FIELD_ENTRIES) {
				printf("MAX_FIELD_ENTRIES exceeded (i=%" PRId32 ", c0=%" PRId32 ", c1=%" PRId32 ".). Ignoring...\n",i,pairs[i].x,pairs[i].y);
				return AMI_ERROR_GENERIC_ERROR;
			}

			qsort(newEntries, numNewEntries, sizeof(s32), cmpS32s);
			int newEndLoc=0;
			int i=0;
			int j=0;
			while (i<numNewEntries) {
				/* How many copies of newEntries[i] are there? */
				j=i+1;
				while ((j<numNewEntries) && (newEntries[j]==newEntries[i])) {
					j++;
				}
				if (mod2) {
					if ((j-i)%2 == 0) {
						/* An even number of them, so omit them altogether. */
						i=j;
					} else {
						/* Copy one of them over. */
						newEntries[newEndLoc++] = newEntries[i];
						i=j;
					}
				} else {
					/* Copy one of them over. */
					newEntries[newEndLoc++] = newEntries[i];
					i=j;
				}
			}
			newBaseStream->write_array(newEntries,newEndLoc);
			releaseField(pairs_y[cdfield++]);
		} else {
			newBaseStream->write_array(buff,fs);
		}
		newOffsets->push_back(newBaseStream->tell());
	}

	delete[] pairs_y; delete[] buff; delete[] pairs_x;
	delete _m_fields; delete _m_baseStream;
	_m_fields=newOffsets; _m_baseStream = newBaseStream;

	return result;
}

AMI_err	_LSList::getSortOnFieldSize(s32 *fields) {
	s32 counts[MAX_LIST_ENTRIES];
	s32 c_idx;
	AMI_err result = AMI_ERROR_NO_ERROR;

	memset(counts,0,sizeof(s32)*MAX_LIST_ENTRIES);
	for (u32 i=0; i<_m_count(); ++i) {
		c_idx=fieldSize(i);
		if (c_idx<MAX_LIST_ENTRIES) {
			++counts[c_idx+1];
		} else {
			printf("getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
			return AMI_ERROR_GENERIC_ERROR;
		}
	}
	for (u32 i=1; i<MAX_LIST_ENTRIES; ++i) {
		counts[i]+=counts[i-1];
	}

	for (u32 i=0; i<_m_count(); ++i) {
		c_idx=fieldSize(i);
		if (c_idx<MAX_LIST_ENTRIES) {
			fields[counts[c_idx]++]=i;
		} else {
			printf("getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
			return AMI_ERROR_GENERIC_ERROR;
		}
	}
	return result;
}

s32 _LSList::getMaxEntry() {
	s32 N = 0;
	s32* item;

	_m_baseStream->seek(0);
	while (_m_baseStream->read_item(&item)!=AMI_ERROR_END_OF_STREAM) {
		if (N < (*item)) {
			N = *item;
		}
	}
	return N;
}

void _LSList::sortEntries(u32 field) {
	TPIE_OS_OFFSET cSize = fieldSize(field);
	s32* entries = new s32[cSize];
	_m_baseStream->seek(_m_fields->at(field));
	_m_baseStream->read_array(entries,&cSize);
	qsort(entries, cSize, sizeof(s32), cmpS32s);
	_m_baseStream->seek(_m_fields->at(field));
	_m_baseStream->write_array(entries,cSize);
	delete[] entries;
}

void _LSList::toLList_t(llist_t* L) {
	ll_init(L,Weight(),numFields());
	L->numFields=numFields();
	for (u32 i=0; i<_m_count(); ++i) {
		_LSListField* cField = fetchField(i);
		for (u32 j=0; j<cField->numEntries();++j) {
			L->data[L->index[i]+j]=cField->at(j);
		}
		L->index[i+1]=L->index[i]+cField->numEntries();
		releaseField(cField);
	}
}

_LookupTable::_LookupTable(s32* offsets, u32 numEntries) {
	_m_entries = new std::vector<s32>;
	_m_temp_entries = new std::vector<s32>;
	data = new s32[offsets[numEntries-1]];
	memset(data,0,offsets[numEntries-1]*sizeof(s32));
	_m_entries->push_back(0);
	_m_temp_entries->push_back(0);
	for (u32 i=0; i<numEntries; ++i) {
		_m_entries->push_back(offsets[i]);
		_m_temp_entries->push_back(offsets[i]);
	}
}

_LookupTable::~_LookupTable() {
	delete _m_entries;
	delete _m_temp_entries;
	delete[] data;
}

_LSListField* _LookupTable::fetchField(u32 field) {
	u32 fs = _m_entries->at(field+1)-_m_entries->at(field);
	s32* entries = new s32[fs];
	memset(entries,0,fs*sizeof(s32));
	for (u32 i=0; i<fs; ++i) {
		entries[i]=data[_m_entries->at(field)+i];
	}
	_LSListField* F = new _LSListField(entries,fs);
	delete[] entries;
	return F;
}

void _LookupTable::releaseField(_LSListField* F) {
	delete F;
}
