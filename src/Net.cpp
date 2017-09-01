#include "Net.h"
#include "Exception.h"

using namespace std;

Net::Net(double tw)
: tilewidth(tw), G((int)(2.0/tw)), gate_set(0), gate_orders(0),
	gate_inverses(0), gate_labels(0)
{
	if (sizeof(int) != 4)
	{
		throw Exception("non 32-bit integers detected - please see Net.h");
		return;
	}
	if (G > 0x7FFE)
	{
		throw Exception("tilewidth too small. (minimum is 6.2 x 10^-5)");
		return;
	}

	memset(label_indices,-1,sizeof(char)*512);
}

Net::Net(const char* file) 
{
	memset(label_indices,-1,sizeof(char)*512);
	load(file);
}


Net::~Net(void)
{
	list<knot*>::iterator ki = knots.begin();
	while (ki != knots.end())
	{
		delete *ki++;
	}

	if (gate_set)
		delete[] gate_set;
	if (gate_orders)
		delete[] gate_orders;
	if (gate_labels)
		delete[] gate_labels;
	if (gate_inverses)
		delete[] gate_inverses;
}


void Net::generate(char* labels, cdouble** basis, int n, int max_length)
{
	// At most we'll need double the number of user supplied gates 
	// (the other ones are inverses)
	
	if (gate_set)
		delete[] gate_set;
	if (gate_orders)
		delete[] gate_orders;
	if (gate_labels)
		delete[] gate_labels;
	if (gate_inverses)
		delete[] gate_inverses;

	gate_set = new Matrix<cdouble>[2*n];
	gate_labels = new char[2*n];
	gate_orders = new int[2*n];
	gate_inverses = new int[2*n];

	N = n;

	// Exact (sort of) equivalence functor
	su2_equiv is_equiv(1e-15);

	extern Matrix<cdouble> Id2;

	// start by adding the supplied gates and checking they're in SU(2)
	for (int i=0;i<N;i++)
	{
		gate_set[i] = Matrix<cdouble>(basis[i],2,2);
		gate_labels[i] = labels[i];

		if (abs(Determinant(gate_set[i])-1.0) > 1e-15)
		{
			char msg[128];
			sprintf(msg, "det(%c) != 1", labels[i]);
			throw Exception(msg);
			return;
		}

		if (!is_equiv(Adjoint(gate_set[i])*gate_set[i],Id2))
		{
			char msg[128];
			sprintf(msg, "%c is not unitary", labels[i]);
			throw Exception(msg);
			return;
		}
	}

	// now add their inverses
	for (int i=0;i<n;i++)
	{
		Matrix<cdouble> GI = Adjoint(gate_set[i]);

		// check it isn't already there
		int j = 0;
		while (j < N && !is_equiv(GI,gate_set[j++]))
			j++;

		// add the matrix to the set if necessary and update the
		// inverses array
		if (j < N)
		{
			gate_inverses[i] = j-1;
			gate_inverses[j-1] = i;
		}
		else
		{
			gate_set[N] = GI;
			gate_labels[N] = gate_labels[i] + 32;	// (lowercase)
			gate_inverses[N] = i;
			gate_inverses[i] = N;
			N++;
		}
	}

	// And match each label to its index in all these arrays
	for (int i=0;i<N;i++)
		label_indices[(int)gate_labels[i]] = i;

	// finally we need to work out the orders. 
	// for our purposes we'll take any order over 50 to be infinite
	for (int i=0;i<N;i++)
	{
		int n = 0;
		Matrix<cdouble> C = gate_set[i];

		while(n++ < 50 && !is_equiv(C,Id2))
			C *= gate_set[i];

		if (n == 50)
			gate_orders[i] = -1;
		else
			gate_orders[i] = n;
	}

	// And fill up the enet
	generate_net(max_length);
}


// This is called from Net::generate above to enumerate all possible sequences
// from the gate set and store them appropriately in the enet.
void Net::generate_net(int max_depth)
{
	int* sequence = new int[max_depth];
	char* _word = new char[max_depth+1];
	char* word = _word+1; _word[0] = '_';
	memset(word,0,max_depth);

	Matrix<cdouble>* products = new Matrix<cdouble>[max_depth];

	memset(sequence,-1,max_depth*sizeof(int));

	int depth = 0;

	while(depth != -1)
	{
		sequence[depth]++;

		if (sequence[depth] < N)	// Append an extra symbol 
		{
			// --- check we're not following something with its inverse
			if (sequence[depth-1] == gate_inverses[sequence[depth]])
				continue;
			// --- check we're not exceeding the order
			if (gate_orders[sequence[depth]] > 0)
			{
				int rpt = 1;
				while (word[depth-rpt] == gate_labels[sequence[depth]]
						&& rpt++ < gate_orders[sequence[depth]])
					;
				if (rpt == gate_orders[sequence[depth]])
					continue;
			}
			// --- 
			
			// --- And here's where we do stuff.
			products[depth] = (depth > 0)
				? products[depth-1]*gate_set[sequence[depth]]
				: gate_set[sequence[depth]];

			word[depth] = gate_labels[sequence[depth]];

			// we're not checking for duplicates here which is silly

			add(products[depth],word,depth+1);

			// ---
			if (depth < max_depth-1)	// go down
			{
				depth++;
			}
		}
		else	// go up 
		{
			word[depth] = 0;
			sequence[depth] = -1;
			depth--;
		}
	}

	delete[] sequence;
	delete[] _word;
	delete[] products;
}


// Map the SU(2) matrix M to a corner of its 'tile'  
Net::icoord Net::coordinate(const Matrix<cdouble>& M, int corner)
{
	icoord MC;
	double mc[4] = {real(M(0,0)),-1.0*imag(M(0,1)),real(M(1,0)),imag(M(1,1))};

	// now turn this real 4-D coordinate into an icoord
	
	// If the corner is -1 (default argument) then we want to return the 
	// nearest corner.
	if (corner < 0)
	{
		MC.t = (short) floor(G*mc[0]+0.5);
		MC.x = (short) floor(G*mc[1]+0.5);
		MC.y = (short) floor(G*mc[2]+0.5);
		MC.z = (short) floor(G*mc[3]+0.5);
	}

	// Otherwise return the specified corner
	else
	{
		if (corner & 1)
			MC.t = (short) ceil(G*mc[0]);
		else
			MC.t = (short) floor(G*mc[0]);

		if (corner & 2)
			MC.x = (short) ceil(G*mc[1]);
		else
			MC.x = (short) floor(G*mc[1]);

		if (corner & 4)
			MC.y = (short) ceil(G*mc[2]);
		else
			MC.y = (short) floor(G*mc[2]);

		if (corner & 8)
			MC.z = (short) ceil(G*mc[3]);
		else
			MC.z = (short) floor(G*mc[3]);
	}

	// We want the first non-zero coordinate to be > 0, so multiply by -I if 
	// necessary.

	if (MC.t < 0)
	{
		MC.t *= -1;
		MC.x *= -1;
		MC.y *= -1;
		MC.z *= -1;
	}
	else if (MC.t == 0)
	{
		if (MC.x < 0)
		{
			MC.x *= -1;
			MC.y *= -1;
			MC.z *= -1;
		}
		else if (MC.x == 0)
		{
			if (MC.y < 0)
			{
				MC.y *= -1;
				MC.z *= -1;
			}
			else if (MC.y == 0)
			{
				if (MC.z < 0)
				{
					MC.z *= -1;
				}
			}
		}
	}

	// MC.t += G: MC.x += G: MC.y += G: MC.z += G:
	return MC;
}

// Basically as above but returns the entire knot
Net::knot* Net::nearest(const Matrix<cdouble>& U)
{
	double dmin = 4.0;
	knot* nn;

	icoord uc = coordinate(U);

	if (su2net[uc].begin() == su2net[uc].end())
	{
		// Gah. Check the 16 surrounding coordinates until we find one
		bool found_one = false;
		for (int i=0;i<16;i++)
		{
			uc = coordinate(U,i);
			if (su2net[uc].begin() != su2net[uc].end())
			{
				found_one = true;
				break;
			}
		}
		if (!found_one)
		{
			throw Exception("unable to find approximation"
					" (try increasing the tile width or sequence length)");
			return 0;
		}
	}

	search_count = 0;
	list<knot*>::iterator ki = su2net[uc].begin();
	while (ki != su2net[uc].end())
	{
		double d = su2::proj_trace_dist(U,(*ki)->M);
		if (d < dmin)
		{
			nn = *ki;
			dmin = d;
		}
		ki++;
		search_count++;
	}

	return nn;
}


// When we add an approximating sequence we need to allocate some memory for it
// (owned by the pointer in knots), and append this pointer to the lists at the
// 16 surrounding coordinates.
//
// There's a problem with the way we're dealing with -1 phases. Imagine we have
// two matrices with coordinates (0,x',y',z') and (-t,x,y,z) with t very small 
// and (x',y',z') ~ (x,y,z). The function Net::coordinate will map the second
// one to (t,-x,-y,-z) which is no-where near the first one.
void Net::add(Matrix<cdouble>& U, char* word, int l)
{
	knot* kn = new knot;
	kn->word.assign(word,l);
	kn->l = l;
	kn->M = U;

	//su2::mat_to_cart4(U,kn->coord);

	knots.push_back(kn);

	for (int c=0;c<16;c++)
	{
		icoord uci = coordinate(U,c);
		su2net[uci].push_back(kn);
	}
}


// Evaluate a string of labels and return the corresponding matrix 
Matrix<cdouble> Net::evaluate(const std::string& s)
{
	int l = s.length();

	Matrix<cdouble> U(2,2);

	if (label_indices[(int)s[0]] == -1)
	{
		throw Exception("unrecognized label in gate sequence");
		return U;
	}

	U = gate_set[label_indices[(int)s[0]]];

	for (int i=1;i<l;i++)
	{
		if (label_indices[(int)s[i]] > N-1)
		{
			throw Exception("unrecognized label in gate sequence");
			return U;
		}
		U *= gate_set[label_indices[(int)s[i]]];
	}

	return U;
}


int Net::load(const char* file)
{
	char header[17];
	char check[] = "Epsilon Net v0.1";
	char msg[128];

	ifstream infile;
	infile.open(file, ios::in | ios::binary);
	if (!infile)
	{
		// let's hope the filename isn't too long
		sprintf(msg,"unable to open file '%s' for reading",file);
		throw Exception(msg);
		return 0;
	}

	// Check that the file actually stores a Net 
	infile.read((char *)header,sizeof(check));
	if (strncmp(header,check,15*sizeof(char)))
	{
		sprintf(msg,"'%s' doesn't look like a Net file",file);
		throw Exception(msg);
		return 0;
	}

	// read in all the structural stuff
	infile.read((char *)&tilewidth,sizeof(tilewidth));
	infile.read((char *)&G,sizeof(G));
	infile.read((char *)&N,sizeof(N));

	gate_set = new Matrix<cdouble>[N];
	gate_labels = new char[N];
	gate_orders = new int[N];
	gate_inverses = new int[N];

	cdouble aa[4];
	for (int i=0;i<N;i++)
	{
		gate_set[i] = Matrix<cdouble>(2,2);
		infile.read((char *)gate_set[i].A,4*sizeof(cdouble));
	}

	infile.read((char *)gate_labels,N*sizeof(char));
	infile.read((char *)gate_orders,N*sizeof(int));
	infile.read((char *)gate_inverses,N*sizeof(int));

	// Match each label to its index in gate_labels etc
	for (int i=0;i<N;i++)
		label_indices[(int)gate_labels[i]] = i;

	// Now read in the actual net
	int total;
	infile.read((char *)&total,sizeof(int));

	char tmp[128];
	double dcoord[4];

	for (int i=0;i<total;i++)
	{
		knot* kn = new knot;
		infile.read((char *)&(kn->l),sizeof(int));
		infile.read((char *)&tmp,kn->l*sizeof(char));
		kn->word.assign(tmp,kn->l);
		infile.read((char *)&(dcoord),4*sizeof(double));
		su2::cart4_to_mat(dcoord,kn->M);

		icoord ic, id;
		for (int c=0;c<16;c++)
		{
			infile.read((char *)&ic,sizeof(icoord));
			su2net[ic].push_back(kn);

			id = coordinate(kn->M,c);
			//cout << ic.t << "," << ic.x << "," << ic.y << "," << ic.z << endl;
			//cout << id.t << "," << id.x << "," << id.y << "," << id.z << endl;
			//cout << "-----\n";
		}
		knots.push_back(kn);
	}

	return 1;
}


void Net::save(const char* file)
{
	char header[] = "Epsilon Net v0.1";
	char msg[128];

	ofstream outfile;
	outfile.open(file, ios::out | ios::binary);

	if (!outfile)
	{
		// let's hope the filename isn't too long
		sprintf(msg,"unable to open file '%s' for writing",file);
		throw Exception(msg);
		return;
	}
	outfile.write(header,sizeof(header));

	outfile.write((char *)&tilewidth,sizeof(tilewidth));
	outfile.write((char *)&G,sizeof(G));
	outfile.write((char *)&N,sizeof(N));

	for (int i=0;i<N;i++)
		outfile.write((char *)gate_set[i].A,4*sizeof(cdouble));

	outfile.write((char *)gate_labels,N*sizeof(char));
	outfile.write((char *)gate_orders,N*sizeof(int));
	outfile.write((char *)gate_inverses,N*sizeof(int));

	// Now dump the entire net. Each matrix is dumped as its 4D euclidean
	// coordinate, and then it is followed up with the 16 coordinates that
	// it is associated with.
	int total = knots.size();
	outfile.write((char *)&total,sizeof(int));

	double dcoord[4];

	list<knot*>::iterator ki = knots.begin();
	while (ki != knots.end())
	{
		// Must write the string out explicity as a bunch of chars
		outfile.write((char *)&(*ki)->l,sizeof(int));
		const char* word = (*ki)->word.c_str();

		//cout << (*ki)->l << " " << word << endl;

		outfile.write((char *)word,(*ki)->l*sizeof(char));

		su2::mat_to_cart4((*ki)->M,dcoord);
		outfile.write((char *)&(dcoord),4*sizeof(double));

		for (int c=0;c<16;c++)
		{
			icoord ic = coordinate((*ki)->M,c);
			outfile.write((char*)&ic,sizeof(icoord));
		}
		ki++;
	}

	outfile.close();
}


string Net::invert(const string& s)
{
	int l = s.length() - 1;

	string ret;

	while (l >= 0)
	{
		if (label_indices[(int)s[l]] == -1)
		{
			throw Exception("unrecognized label in gate sequence");
			return "";
		}

		ret += gate_labels[gate_inverses[label_indices[(int)s[l--]]]];
	}

	return ret;
}


ostream& operator<<(ostream& out, const Net& N)
{
	out << "Net of tilewidth " << N.tilewidth << std::endl;
	out << "gate set { ";
	for (int i=0;i<N.N;i++)
	{
		out << N.gate_labels[i] << " (order " << N.gate_orders[i];
		if (N.gate_orders[i] == 2)
			out << ")";
		else
			out << ", inverse " << N.gate_labels[N.gate_inverses[i]] << ")";

		if (i != N.N-1)
			out << ", ";
	}
	out << " }";

	//out << N.knots.size() << " net points";

	return out;
}


Net::knot Net::solovay_kitaev(const Matrix<cdouble>& U, int depth)
{
	if (depth == 0)
	{
		return *nearest(U);
	}

	knot ku = solovay_kitaev(U,depth-1);

	Matrix<cdouble> V(2,2), W(2,2);
	su2::group_factor(U*Adjoint(ku.M),V,W);

	knot kv = solovay_kitaev(V,depth-1);
	knot kw = solovay_kitaev(W,depth-1);

	// Construct the returning knot
	knot ret =
	{
		2*(kv.l) + 2*(kw.l) + ku.l,	
		kv.word + kw.word + invert(kv.word) + invert(kw.word) + ku.word,
		kv.M * kw.M * Adjoint(kv.M) * Adjoint(kw.M) * ku.M
	};

	return ret;
}
