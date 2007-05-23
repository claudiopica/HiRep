template<class INDEX, class VALUE> class sparselinear : public orderedlist<INDEX,VALUE>
{
public:
	int size;
	
	sparselinear() : orderedlist<INDEX,VALUE>() { size = 0; }
	sparselinear(int n) : orderedlist<INDEX,VALUE>() { size = n; }
	sparselinear(const sparselinear<INDEX,VALUE>& mat) : orderedlist<INDEX,VALUE>(mat) { size = mat.size; }
	
	virtual sparselinear<INDEX,VALUE>& operator=(const sparselinear<INDEX,VALUE>& b)
	{
		for(KEYTYPE i = 0; i < orderedlist<INDEX,VALUE>::length; i++)
			delete orderedlist<INDEX,VALUE>::data[i];
		if(orderedlist<INDEX,VALUE>::length != 0) delete[] orderedlist<INDEX,VALUE>::data;
		size = b.size;
		orderedlist<INDEX,VALUE>::length = b.length;
		orderedlist<INDEX,VALUE>::data = new typename orderedlist<INDEX,VALUE>::ELEMENT*[orderedlist<INDEX,VALUE>::length];
		for(KEYTYPE i = 0; i < orderedlist<INDEX,VALUE>::length; i++)
		{
			orderedlist<INDEX,VALUE>::data[i] = new typename orderedlist<INDEX,VALUE>::ELEMENT;
			orderedlist<INDEX,VALUE>::data[i]->index = b.data[i]->index;
			orderedlist<INDEX,VALUE>::data[i]->value = b.data[i]->value;
		}
		return *this;
	}
	
	virtual void minus()
	{
		for(int i = 0; i < orderedlist<INDEX,VALUE>::length; i++)
			(orderedlist<INDEX,VALUE>::data)[i]->value.minus();
	}
	
	virtual void add(sparselinear<INDEX,VALUE>& l)
	{
		if(size != l.size)
		{
			cerr << "WRONG LINEAR SUM.\n";
			return;
		}
		orderedlist<INDEX,VALUE>::add(l);
	}
	virtual void sub(const sparselinear<INDEX,VALUE>& l)
	{
		if(size != l.size)
		{
			cerr << "WRONG LINEAR SUM.\n";
			return;
		}
		minus();
		orderedlist<INDEX,VALUE>::add(l);
		minus();
	}
	
	friend void herm(VALUE& ret, const sparselinear<INDEX,VALUE>& l, const sparselinear<INDEX,VALUE>& k)
	{
		if(k.size != l.size)
		{
			cerr << "WRONG LINEAR HERM.\n";
			return;
		}
		ret.clear();
		for(int i = 0; i < l.length; i++)
			for(int j = 0; j < k.length; j++)
				if(l[i].index == k[j].index)
				{
					VALUE tmp;
					tmp = l[i].value;
					tmp.conjugate();
					ret.add_mult(tmp, k[j].value);
				}
	}
};


struct mindex
{
	int row;
	int col;
	int size;
	
	mindex() { row = 0; col = 0; size = 0;}
	mindex(int r, int c, int s)
	{
		if(s < 1)
		{
			cerr << "mindex: WRONG SIZE.\n";
			size = 0;
		}
		else
			size = s;
		if(r < 0 || r >= size)
		{
			cerr << "mindex: OUT OF RANGE (row = " << r << ", size = " << size << ").\n";
			row = 0;
		}
		else
			row = r;
		if(c < 0 || c >= size)
		{
			cerr << "mindex: OUT OF RANGE (col = " << c << ", size = " << size << ").\n";
			col = 0;
		}
		else
			col = c;
	}
	
	mindex& operator=(const mindex& index)
	{
		row = index.row;
		col = index.col;
		size = index.size;
		return *this;
	}

	int key() const { return row*size+col; }
	friend bool operator<(const mindex& a, const mindex& b) { return a.key() < b.key(); }
	friend bool operator>(const mindex& a, const mindex& b) { return a.key() > b.key(); }
	friend bool operator<=(const mindex& a, const mindex& b) { return a.key() <= b.key(); }
	friend bool operator>=(const mindex& a, const mindex& b) { return a.key() >= b.key(); }
	friend bool operator==(const mindex& a, const mindex& b) { return a.key() == b.key() && a.size == b.size; }
	friend bool operator!=(const mindex& a, const mindex& b) { return a.key() != b.key() || a.size != b.size; }
	
	friend ostream& operator<<(ostream& os, const mindex& a)
	{
		return os << _MINDEX_;
	}
};

template<class VALUE> class sparsematrix : public sparselinear<mindex,VALUE>
{
public:
	sparsematrix() : sparselinear<mindex,VALUE>() {}
	sparsematrix(int n) : sparselinear<mindex,VALUE>(n) {}
	sparsematrix(const sparsematrix<VALUE>& mat) : sparselinear<mindex,VALUE>(mat) {}

	const VALUE& get(int row, int col) const
	{
		return orderedlist<mindex,VALUE>::get(mindex(row,col,sparselinear<mindex,VALUE>::size));
	}
	void set(int row, int col, const VALUE& value)
	{
		orderedlist<mindex,VALUE>::set(mindex(row,col,sparselinear<mindex,VALUE>::size),value);
	}
	void add(int row, int col, const VALUE& value)
	{
		orderedlist<mindex,VALUE>::add(mindex(row,col,sparselinear<mindex,VALUE>::size),value);
	}

	virtual string getstring() const
	{
		ostringstream RET;
		for(int i = 0; i < sparselinear<mindex,VALUE>::size; i++)
		{
			RET << "[ ";
			for(int j = 0; j < sparselinear<mindex,VALUE>::size; j++)
				RET << get(i,j) << " ";
			RET << "]\n";
		}
		RET << "\n";
		return RET.str();
	}
	
	virtual void adjoint()
	{
		int tmp;
		for(int i = 0; i < this->length; i++)
		{
			tmp = this->data[i]->index.row;
			this->data[i]->index.row = this->data[i]->index.col;
			this->data[i]->index.col = tmp;
			this->data[i]->value.conjugate();
		}
		this->sort();
	}
	virtual void transpose()
	{
		int tmp;
		for(int i = 0; i < this->length; i++)
		{
			tmp = this->data[i]->index.row;
			this->data[i]->index.row = this->data[i]->index.col;
			this->data[i]->index.col = tmp;
		}
		this->sort();
	}
	
	virtual void add(sparsematrix<VALUE>& m)
	{
		sparselinear<mindex,VALUE>::add(m);
	}
	virtual void sub(sparsematrix<VALUE>& m)
	{
		sparselinear<mindex,VALUE>::sub(m);
	}
	void mult(const sparsematrix<VALUE>& p, const sparsematrix<VALUE>& q)
	{
		if(p.size != q.size)
		{
			cerr << "WRONG MATRIX MULT.\n";
			return;
		}
		sparselinear<mindex,VALUE>::clear();
		sparselinear<mindex,VALUE>::size = p.size;
		for(int i = 0; i < p.length; i++)
			for(int j = 0; j < q.length; j++)
				if(p[i].index.col == q[j].index.row)
				{
					VALUE value;
					value.mult(p[i].value, q[j].value);
					add(p[i].index.row, q[j].index.col, value);
				}
	}
	void add_mult(const sparsematrix<VALUE>& p, const sparsematrix<VALUE>& q)
	{
		if(p.size != q.size || sparselinear<mindex,VALUE>::size != p.size)
		{
			cerr << "WRONG MATRIX MULT.\n";
			return;
		}
		for(int i = 0; i < p.length; i++)
			for(int j = 0; j < q.length; j++)
				if(p[i].index.col == q[j].index.row)
				{
					VALUE value;
					value.mult(p[i].value, q[j].value);
					add(p[i].index.row, q[j].index.col, value);
				}
	}
	friend void trace(VALUE& ret, const sparsematrix<VALUE>& l)
	{
		ret.clear();
		for(int i = 0; i < l.length; i++)
			if(l[i].index.row == l[i].index.col)
				ret.add(l[i].value);
	}
};


struct vindex
{
	int row;
	int size;
	
	vindex() { row = 0; size = 0;}
	vindex(int r, int s)
	{
		if(s < 1)
		{
			cerr << "mindex: WRONG SIZE.\n";
			size = 0;
		}
		else
			size = s;
		if(r < 0 || r >= size)
		{
			cerr << "mindex: OUT OF RANGE (row = " << r << ", size = " << size << ").\n";
			row = 0;
		}
		else
			row = r;
	}
	
	vindex& operator=(const vindex& index)
	{
		row = index.row;
		size = index.size;
		return *this;
	}

	friend bool operator<(const vindex& a, const vindex& b) { return a.row < b.row; }
	friend bool operator>(const vindex& a, const vindex& b) { return a.row > b.row; }
	friend bool operator<=(const vindex& a, const vindex& b) { return a.row <= b.row; }
	friend bool operator>=(const vindex& a, const vindex& b) { return a.row >= b.row; }
	friend bool operator==(const vindex& a, const vindex& b) { return a.row == b.row && a.size == b.size; }
	friend bool operator!=(const vindex& a, const vindex& b) { return a.row != b.row || a.size != b.size; }
	
	friend ostream& operator<<(ostream& os, const vindex& a)
	{
		return os << _VINDEX_;
	}
};

template<class VALUE> class sparsevector : public sparselinear<vindex,VALUE>
{
public:
	sparsevector() : sparselinear<vindex,VALUE>() {}
	sparsevector(int n) : sparselinear<vindex,VALUE>(n) {}
	sparsevector(const sparsevector<VALUE>& vec) : sparselinear<vindex,VALUE>(vec) {}
/*	
	sparsevector<VALUE>& operator=(const sparsevector<VALUE>& b)
	{
		for(KEYTYPE i = 0; i < orderedlist<vindex,VALUE>::length; i++)
			delete orderedlist<vindex,VALUE>::data[i];
		if(orderedlist<vindex,VALUE>::length != 0) delete[] orderedlist<vindex,VALUE>::data;
		size = b.size;
		orderedlist<vindex,VALUE>::length = b.length;
		orderedlist<vindex,VALUE>::data = new typename orderedlist<vindex,VALUE>::ELEMENT*[orderedlist<vindex,VALUE>::length];
		for(KEYTYPE i = 0; i < orderedlist<int,VALUE>::length; i++)
		{
			orderedlist<vindex,VALUE>::data[i] = new typename orderedlist<vindex,VALUE>::ELEMENT;
			orderedlist<vindex,VALUE>::data[i]->index = b.data[i]->index;
			orderedlist<vindex,VALUE>::data[i]->value = b.data[i]->value;
		}
		return *this;
	}
*/
	const VALUE& get(int row) const
	{
		return orderedlist<vindex,VALUE>::get(vindex(row,sparselinear<vindex,VALUE>::size));
	}
	void set(int row, const VALUE& value)
	{
		orderedlist<vindex,VALUE>::set(vindex(row,sparselinear<vindex,VALUE>::size),value);
	}
	void add(int row, const VALUE& value)
	{
		orderedlist<vindex,VALUE>::add(vindex(row,sparselinear<vindex,VALUE>::size),value);
	}

	virtual string getstring() const
	{
		ostringstream RET;
		RET << "[ ";
		for(int i = 0; i < sparselinear<vindex,VALUE>::size; i++)
			RET << get(i) << " ";
		RET << "]\n";
		return RET.str();
	}

	virtual void add(sparsevector<VALUE>& v)
	{
		sparselinear<vindex,VALUE>::add(v);
	}
	virtual void sub(sparsevector<VALUE>& v)
	{
		sparselinear<vindex,VALUE>::sub(v);
	}
	void dot(const sparsematrix<VALUE>& p, const sparsevector<VALUE>& v)
	{
		if(p.size != v.size)
		{
			cerr << "WRONG VECTOR DOT.\n";
			return;
		}
		sparselinear<vindex,VALUE>::clear();
		sparselinear<vindex,VALUE>::size = v.size;
		for(int i = 0; i < p.length; i++)
			for(int j = 0; j < v.length; j++)
				if(p[i].index.col == v[j].index.row)
				{
					VALUE value;
					value.mult(p[i].value, v[j].value);
					add(p[i].index.row, value);
				}
	}
	void add_dot(const sparsematrix<VALUE>& p, const sparsevector<VALUE>& v)
	{
		if(sparselinear<vindex,VALUE>::size != p.size || sparselinear<vindex,VALUE>::size != v.size)
		{
			cerr << "WRONG VECTOR DOT.\n";
			return;
		}
		for(int i = 0; i < p.length; i++)
			for(int j = 0; j < v.length; j++)
				if(p[i].index.col == v[j].index.row)
				{
					VALUE value;
					value.mult(p[i].value, v[j].value);
					add(p[i].index.row, value);
				}
	}

};
