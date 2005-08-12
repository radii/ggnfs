/**************************************************************
 * llist_tpie.h
 * Copyright 2005, Anton Korobeynikov, Chris Monico.
 **************************************************************/

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

#ifndef __LList_TPIE_H__
#define __LList_TPIE_H__

#include "ggnfs.h"
#include "app_config.h"

#include <iostream>

#include <utility>
#include <list>
#include <vector>

#include <assert.h>

#include <ami_stream.h>
#include <ami_block.h>
#include <ami_coll.h>
#include <ami_stack.h>
#include <ami_cache.h>

#define BLOCK_SIZE 50*1024
#define MAX_FIELD_ENTRIES 2048
#define MAX_LIST_ENTRIES 10

typedef TPIE_OS_SIZE_T entry_count_t;

struct _LBList_entry_info {
	entry_count_t entries;
};

int    cmpS32s(const void *a, const void *b);
int    cmpU32s(const void *a, const void *b);
int    cmp_lpair_t(const void *a, const void *b);

template<class BTECOLL> class _LBList;

template<class BTECOLL = BTE_COLLECTION>
class _LBListField: public AMI_block<s32, _LBList_entry_info, BTECOLL> {
public:

	using AMI_block<s32, _LBList_entry_info, BTECOLL>::info;
	using AMI_block<s32, _LBList_entry_info, BTECOLL>::el;
	using AMI_block<s32, _LBList_entry_info, BTECOLL>::lk;
	using AMI_block<s32, _LBList_entry_info, BTECOLL>::dirty;
	using AMI_block<s32, _LBList_entry_info, BTECOLL>::bid;

	_LBListField(AMI_collection_single<BTECOLL>* pcoll, AMI_bid ID,	entry_count_t maxEntries=0);
	_LBListField(const _LBListField<BTECOLL>& other);
	_LBListField<BTECOLL>& operator=(const _LBListField<BTECOLL>& other);

protected:
	entry_count_t _maxEntries;
private:
	//  No private members.
};

template<class BTECOLL = BTE_COLLECTION>
class _LBList {
public:
	_LBList(size_t lbf = 1, persistence per = PERSIST_DELETE);
	_LBList(llist_t* L, size_t lbf = 1, persistence per = PERSIST_DELETE);
	~_LBList();

	static const size_t	getLBFbyEntries (size_t maxEntries) { return ((maxEntries*sizeof(s32)+sizeof(entry_count_t))/TPIE_OS_BLOCKSIZE()+1); };
	void			clearList();
	inline unsigned int	numFields() const { return _count; };
	_LBListField<BTECOLL>*	fetchField(unsigned int index);
	void			releaseField(_LBListField<BTECOLL>* field);

	AMI_err			appendField(s32 *entries, s32 numEntries);
	AMI_err			appendField(const _LBListField<BTECOLL>& newEntry);
	AMI_err			appendToField(s32 field, s32 *entries, s32 numEntries);
	AMI_err			deleteFields(s32 *fields, s32 numFields);
	int				numCommonEntries(s32 field0, s32 field1);
	AMI_err			catFields(lpair_t *pairs, s32 numPairs, int mod2);
	AMI_err			getSortOnFieldSize(s32 *fields);
	llist_t*		toLList_t();
	s32				getMaxEntry();
	s32				fieldSize(u32 field);
	s32				Weight();
	void			sortEntries(u32 field);

protected:
	AMI_bid* 						_entries;
	unsigned int 					_count;
	entry_count_t 					_maxEntries;
	AMI_collection_single<BTECOLL>*	_collection;

	class remove_entry { 
	public:
		void operator()(_LBListField<BTECOLL>* p) { delete p; }
	};

	typedef AMI_CACHE_MANAGER<_LBListField<BTECOLL>*, remove_entry> entry_cache_t;

	entry_cache_t* _entry_cache;

private:
	void reallocateEntries();
};

class _LSListField {
public:
	_LSListField(s32 *entries=NULL, s32 numEntries=0);
	virtual	~_LSListField();

	inline s32& 	operator[] (u32 entry) const { return _m_entries->at(entry); };
	inline s32& 	at (u32 entry) const { return _m_entries->at(entry); }
	inline u32		numEntries() const { return _m_entries->size(); };
	void			appendEntry(s32 entry);

protected:
	std::vector<s32>* _m_entries;
};

class _LSList {

public:
	_LSList(persistence per = PERSIST_DELETE);
	_LSList(llist_t* L, persistence per = PERSIST_DELETE);
	~_LSList();

	inline u32			numFields() const { return _m_fields->size()-1; };
	_LSListField*		fetchField(u32 index);
	void				releaseField(_LSListField* field);

	AMI_err				appendField(s32 *entries, s32 numEntries);
	AMI_err				appendField(const _LSListField& newField);
	AMI_err				deleteFields(s32 *fields, s32 numFields);
	int					numCommonEntries(s32 field0, s32 field1);
	AMI_err				catFields(lpair_t *pairs, s32 numPairs, bool mod2);
	AMI_err				getSortOnFieldSize(s32 *fields);
	void				toLList_t(llist_t* L);
	s32					getMaxEntry();
	inline s32			fieldSize(u32 field) const { return (_m_fields->at(field+1)-_m_fields->at(field)); };
	inline s32			Weight() const { return _m_fields->at(_m_fields->size()-1); };
	void				sortEntries(u32 field);

protected:
	std::vector<off_t>*			_m_fields;
	AMI_stream<s32>* 			_m_baseStream;
	persistence 				_m_per;
	inline u32					_m_count() const { return _m_fields->size()-1; };
};

class _LookupTable {
public:
	_LookupTable(s32* offsets, u32 numEntries);
	virtual ~_LookupTable();

	inline u32			numFields() const { return _m_entries->size()-1; };
	_LSListField*		fetchField(u32 field);
	void				releaseField(_LSListField* F); 
	inline void			appendEntry(u32 field, s32 entry) {	data[(_m_temp_entries->at(field))++]=entry; };

protected:
	// we don't care about memory here. Just about speed & safety!
	std::vector<s32>*	_m_entries;
	std::vector<s32>*	_m_temp_entries;
	s32* 				data;
};

#endif
