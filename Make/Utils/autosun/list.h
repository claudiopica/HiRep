typedef long long int KEYTYPE;

template<class INDEX, class VALUE> class orderedlist
{
protected:
	struct ELEMENT
	{
		INDEX index;
		VALUE value;
	};
	
	ELEMENT** data;
	
	void sort(int l, int r)
	{
		KEYTYPE left = l, right = r;
		KEYTYPE flag = left+(right-left+1)/2;
		ELEMENT* tmp;

		while(left != right)
		{
			if(data[left]->index <= data[flag]->index && left != flag)
				left++;
			else if(data[flag]->index <= data[right]->index && right != flag)
				right--;
			else if(left != flag && right != flag)
			{
				memcpy((void*)(&tmp), (void*)(data+left), sizeof(ELEMENT*));
				memcpy((void*)(data+left), (void*)(data+right), sizeof(ELEMENT*));
				memcpy((void*)(data+right), (void*)(&tmp), sizeof(ELEMENT*));
			}
			else if(left == flag && right != flag)
			{
				memcpy((void*)(&tmp), (void*)(data+left), sizeof(ELEMENT*));
				memcpy((void*)(data+left), (void*)(data+right), sizeof(ELEMENT*));
				memcpy((void*)(data+right), (void*)(&tmp), sizeof(ELEMENT*));
				flag = right;
				left++;
			}
			else if(left != flag && right == flag)
			{
				memcpy((void*)(&tmp), (void*)(data+left), sizeof(ELEMENT*));
				memcpy((void*)(data+left), (void*)(data+right), sizeof(ELEMENT*));
				memcpy((void*)(data+right), (void*)(&tmp), sizeof(ELEMENT*));
				flag = left;
				right--;
			}
		}

		if(flag-1 > l)
			sort(l, flag-1);
		if(flag+1 < r)
			sort(flag+1, r);
	}


public:
	KEYTYPE length;
	
	orderedlist()
	{
		length = 0;
	}
	orderedlist(const INDEX& index, const VALUE& value)
	{
		length = 1;
		data = new ELEMENT*[1];
		data[0] = new ELEMENT;
		data[0]->index = index;
		data[0]->value = value;
	}
	orderedlist(const orderedlist<INDEX,VALUE>& list)
	{
		length = list.length;
		data = new ELEMENT*[length];
		for(KEYTYPE i = 0; i < length; i++)
		{
			data[i] = new ELEMENT;
			data[i]->index = list.data[i]->index;
			data[i]->value = list.data[i]->value;
		}
		cerr << "COPY CONSTRUCTOR.\n";
	}
	virtual ~orderedlist() = 0;
	
	ELEMENT& operator[](KEYTYPE i) const { return *(data[i]); }
	virtual orderedlist<INDEX,VALUE>& operator=(const orderedlist<INDEX,VALUE>& list)
	{
		for(KEYTYPE i = 0; i < length; i++)
			delete data[i];
		if(length != 0) delete[] data;
		length = list.length;
		data = new ELEMENT*[length];
		for(KEYTYPE i = 0; i < length; i++)
		{
			data[i] = new ELEMENT;
			data[i]->index = list.data[i]->index;
			data[i]->value = list.data[i]->value;
		}
		return *this;
	}
	bool operator==(const orderedlist<INDEX,VALUE>& list) const
	{
		if(length != list.length) return false;
		for(KEYTYPE i = 0; i < length; i++)
		{
			if(data[i]->index != list.data[i]->index) return false;
			if(data[i]->value != list.data[i]->value) return false;
		}
		return true;
	}
	bool operator!=(const orderedlist<INDEX,VALUE>& list) const
	{
		if(length != list.length) return true;
		for(KEYTYPE i = 0; i < length; i++)
		{
			if(data[i]->index != list.data[i]->index) return true;
			if(data[i]->value != list.data[i]->value) return true;
		}
		return false;
	}
	virtual void add(const orderedlist<INDEX,VALUE>& list)
	{
		for(int i = 0; i < list.length; i++)
			add(list[i].index,list[i].value);
	}
	virtual void scale(const VALUE& factor)
	{
		for(int i = 0; i < length; i++)
			data[i]->value *= factor;
	}
	
	void sort()
	{
		sort(0, length-1);
	}
	
	void clear()
	{
		for(KEYTYPE i = 0; i < length; i++)
			delete data[i];
		if(length != 0) delete[] data;
		length = 0;
	}
	virtual VALUE& getzero() const = 0;
	const VALUE& get(const INDEX& index) const
	{
		if(length == 0)
			return getzero();
		else if(data[length-1]->index < index)
			return getzero();
		else
		{
			KEYTYPE pos;
			KEYTYPE start = 0;
			KEYTYPE end = length-1;
			while(1)
			{
				pos = start + (end-start)/2;
				if(index == data[pos]->index || start == end) break;
				else if(index <= data[pos]->index) end = pos;
				else start = pos+1;
			}
			
			if(index == data[pos]->index)
				return data[pos]->value;
			else
				return getzero();
		}
	}
	void set(const INDEX& index, const VALUE& value)
	{
		if(value == getzero())
		{
			remove(index);
			return;
		}
		
		if(length == 0)
		{
			length = 1;
			data = new ELEMENT*[1];
			data[0] = new ELEMENT;
			data[0]->index = index;
			data[0]->value = value;
		}
		else if(data[length-1]->index < index)
		{
			ELEMENT** olddata = data;
			length++;
			data = new ELEMENT*[length];
			memcpy((void*)data, (void*)olddata, sizeof(ELEMENT*)*(length-1));
			delete[] olddata;
			data[length-1] = new ELEMENT;
			data[length-1]->index = index;
			data[length-1]->value = value;
		}
		else
		{
			KEYTYPE pos;
			KEYTYPE start = 0;
			KEYTYPE end = length-1;
			while(1)
			{
				pos = start + (end-start)/2;
				if(index == data[pos]->index || start == end) break;
				else if(index <= data[pos]->index) end = pos;
				else start = pos+1;
			}
			
			if(index == data[pos]->index)
			{
				data[pos]->value = value;
			}
			else if(pos == 0)
			{
				ELEMENT** olddata = data;
				length++;
				data = new ELEMENT*[length];
				memcpy((void*)(data+1), (void*)olddata, sizeof(ELEMENT*)*(length-1));
				delete[] olddata;
				data[pos] = new ELEMENT;
				data[pos]->index = index;
				data[pos]->value = value;
			}
			else
			{
				ELEMENT** olddata = data;
				length++;
				data = new ELEMENT*[length];
				memcpy((void*)data, (void*)olddata, sizeof(ELEMENT*)*pos);
				memcpy((void*)(data+pos+1), (void*)(olddata+pos), sizeof(ELEMENT*)*(length-pos-1));
				delete[] olddata;
				data[pos] = new ELEMENT;
				data[pos]->index = index;
				data[pos]->value = value;
			}
		}
	}
	void add(const INDEX& index, const VALUE& value)
	{
		if(value == getzero())
			return;
		
		if(length == 0)
		{
			length = 1;
			data = new ELEMENT*[1];
			data[0] = new ELEMENT;
			data[0]->index = index;
			data[0]->value = value;
		}
		else if(data[length-1]->index < index)
		{
			ELEMENT** olddata = data;
			length++;
			data = new ELEMENT*[length];
			memcpy((void*)(data), (void*)olddata, sizeof(ELEMENT*)*(length-1));
			delete[] olddata;
			data[length-1] = new ELEMENT;
			data[length-1]->index = index;
			data[length-1]->value = value;
		}
		else
		{
			KEYTYPE pos;
			KEYTYPE start = 0;
			KEYTYPE end = length-1;
			while(1)
			{
				pos = start + (end-start)/2;
				if(index == data[pos]->index || start == end) break;
				else if(index <= data[pos]->index) end = pos;
				else start = pos+1;
			}
			
			if(index == data[pos]->index)
			{
				data[pos]->value += value;
				if(data[pos]->value == getzero())
				{
					ELEMENT** olddata = data;
					length--;
					data = new ELEMENT*[length];
					memcpy((void*)data, (void*)olddata, sizeof(ELEMENT*)*pos);
					memcpy((void*)(data+pos), (void*)(olddata+pos+1), sizeof(ELEMENT*)*(length-pos));
					delete[] olddata;
				}
			}
			else if(pos == 0)
			{
				ELEMENT** olddata = data;
				length++;
				data = new ELEMENT*[length];
				memcpy((void*)(data+1), (void*)olddata, sizeof(ELEMENT*)*(length-1));
				delete[] olddata;
				data[pos] = new ELEMENT;
				data[pos]->index = index;
				data[pos]->value = value;
			}
			else
			{
				ELEMENT** olddata = data;
				length++;
				data = new ELEMENT*[length];
				memcpy((void*)data, (void*)olddata, sizeof(ELEMENT*)*pos);
				memcpy((void*)(data+pos+1), (void*)(olddata+pos), sizeof(ELEMENT*)*(length-pos-1));
				delete[] olddata;
				data[pos] = new ELEMENT;
				data[pos]->index = index;
				data[pos]->value = value;
			}
		}
	}
	void remove(const INDEX& index)
	{
		if(length == 0)
			return;
		else if(data[length-1]->index < index)
			return;
		else
		{
			KEYTYPE pos;
			KEYTYPE start = 0;
			KEYTYPE end = length-1;
			while(1)
			{
				pos = start + (end-start)/2;
				if(index == data[pos]->index || start == end) break;
				else if(index <= data[pos]->index) end = pos;
				else start = pos+1;
			}
			
			if(index == data[pos]->index)
			{
				ELEMENT** olddata = data;
				length--;
				data = new ELEMENT*[length];
				memcpy((void*)data, (void*)olddata, sizeof(ELEMENT*)*pos);
				memcpy((void*)(data+pos), (void*)(olddata+pos+1), sizeof(ELEMENT*)*(length-pos));
				delete[] olddata;
			}
		}
	}

	virtual string getstring() const = 0;
	friend ostream& operator<<(ostream& os, const orderedlist<INDEX,VALUE>& list)
	{
		os << list.getstring();
		return os;
	}
	void printaslist() const
	{
		cout << this << "\t(length = " << length << ")\n";
		for(int i = 0; i < length; i++)
			cout << data[i]->index << "\t" << data[i]->value << "\n";
	}
	void print() const
	{
		cout << this << "\t(length = " << length << ")\n";
		cout << getstring();
	}
};

template<class INDEX, class VALUE> orderedlist<INDEX,VALUE>::~orderedlist<INDEX,VALUE>()
{
	for(KEYTYPE i = 0; i < length; i++)
		delete data[i];
	if(length != 0) delete[] data;
}
