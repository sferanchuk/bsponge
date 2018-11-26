
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <string>
#include <set>
#include <map>
#include <vector>

using namespace std;

enum { clSingle = 1, clComplete, clUpgma };

typedef pair<int,int> IPair;

const double big = 1e20;

void hierarchicalClustering( vector< vector<float> >& pairwiseDistances, double threshold, vector< set<int> >& clindex, int linkage = clComplete ) 
{
	double minDistance;
	//vector< vector< int > > clindex;
	vector< int > clusters;
	int clusterA, clusterB;
	clindex.resize( pairwiseDistances.size() );
	//rclusters.resize( pairwiseDistances.size() );
	for ( int i = 0; i < pairwiseDistances.size(); i++ ) 
	{
		clindex[i].insert( i );
		clusters.push_back( i );
		//rclusters[ i ].insert( i );
	}
	int rcnt = pairwiseDistances.size() + 1;
	while ( true ) 
	{
		if ( clindex.size() == 1 ) 
		{
			break;	
		}
		double minDistance = big;
		IPair candidates;	

		for ( int i = 0; i < clindex.size(); i++ ) 
		{
			clusterA = i;
			for ( int j = i + 1; j < clindex.size(); j++ ) 
			{
				clusterB = j;
				double distance = pairwiseDistances[ clusterA ][ clusterB ];
				if ( distance < minDistance) 
				{
					minDistance = distance;
					candidates = IPair( clusterA, clusterB );
				}
			}
		}
		if ( minDistance > threshold )
		{
			break;
		}
		int sizeA = clindex[ candidates.first ].size();
		int sizeB = clindex[ candidates.second ].size();
		int size = sizeA + sizeB;
		double distance = minDistance;
		set<int> nindex = clindex[ candidates.first ];
		nindex.insert( clindex[ candidates.second ].begin(), clindex[ candidates.second ].end() );
		clindex[ candidates.first ] = nindex;
		clindex.erase( clindex.begin() + candidates.second );
		for ( int i = 0; i < clindex.size(); i++ ) 
		{
			int ii = i;
			if ( ii >= candidates.second ) ii++;
			if ( i == candidates.first ) 
			{
				distance = 0;
			} 
			else 
			{
				double distanceA = pairwiseDistances[ candidates.first ][ ii ];
				double distanceB = pairwiseDistances[ candidates.second ][ ii ];
				switch ( linkage ) {
					case 1:
						distance = fmin( distanceA, distanceB ); 
						break;
					case 2:
						distance = fmax( distanceA, distanceB ); 
						break;
					case 3:
						distance = ( ( distanceA * sizeA ) + ( distanceB * sizeB ) ) / size;
						break;
				}
			}
			pairwiseDistances[ candidates.first ][ ii ] = distance;
			pairwiseDistances[ ii ][ candidates.first ] = distance;
		}
		for ( int i = 0; i < pairwiseDistances.size(); i++ ) 
		{
			pairwiseDistances[ i ].erase( pairwiseDistances[ i ].begin() + candidates.second );
		}
		// Remove row of second cluster from pairwise distance matrix.
		pairwiseDistances.erase( pairwiseDistances.begin() + candidates.second );
	}
}

bool CalcDistances( vector<string>& align, vector< vector<float> >& distances )
{
	distances.resize( align.size() );
	for ( int sc = 0; sc < align.size(); sc++ )
	{
		distances[sc].resize( align.size(), 1. );
		for ( int sc1 = 0; sc1 <= sc; sc1++ )
		{
			if ( align[sc1].size() != align[sc].size() )
			{
				fprintf( stderr, "sequences are not aligned\n" );
				return false;
			}
			int smuts = 0;
			int pg[2] = {0};
			int ml[2] = {0};
			for ( int pos = 0; pos < align[ sc ].size(); pos++ )
			{
				if ( align[sc][pos] == '-' && align[sc1][pos] =='-' ) continue;
				if ( align[sc][pos] == '-' )
				{
					pg[1] = 0;
					if ( pg[0] ) continue;
					else 
					{
						smuts++;
						pg[0] = 1;
						ml[1]++;
					}
				}
				else if ( align[sc1][pos] == '-' )
				{
					pg[0] = 0;
					if ( pg[1] ) continue;
					else 
					{
						smuts++;
						pg[1] = 1;
						ml[0]++;
					}
				}
				else
				{
					ml[0]++;
					ml[1]++;
					pg[0] = 0;
					pg[1] = 0;
					if ( toupper( align[sc][pos] ) != toupper( align[ sc1 ][ pos ] ) ) smuts++;
				}
			}
			int size = min( ml[0], ml[1] );
			double dist = double( smuts ) / size;
			distances[ sc ][ sc1 ] = dist;
		}
	}
	for ( int sc = 0; sc < align.size(); sc++ )
	{
		for ( int sc1 = 0; sc1 < sc; sc1++ )
		{
			distances[ sc1 ][ sc ] = distances[ sc ][ sc1 ];
		}
	}
	return true;
}

const size_t alphabets = 26;

int min3( int a, int b, int c )
{
    if (a <= b && a <= c)
        return a;
    else if (b <= a && b <= c)
        return b;
    else
        return c;
}

int NWAlign(const string &a, const string &b, int alpha_gap, int alpha[alphabets][alphabets], string &a_aligned, string &b_aligned )
{
    size_t n = a.size();
    size_t m = b.size();

    vector<vector<int> > A(n + 1, vector<int>(m + 1));

    for (size_t i = 0; i <= m; ++i)
        A[0][i] = alpha_gap * i;
    for (size_t i = 0; i <= n; ++i)
        A[i][0] = alpha_gap * i;

    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];
            A[i][j] = min3( A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'],
                            A[i-1][j] + alpha_gap,
                            A[i][j-1] + alpha_gap );
        }
    }

    // print2DVector(A);

    a_aligned = "";
    b_aligned = "";
    size_t j = m;
    size_t i = n;
    for (; i >= 1 && j >= 1; )
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j] == A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'])
        {
            /*
             * I think prepending chars this way to a std::string is very inefficient.
             * Is there any better way of doing this without using C-style strings?
             */
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
			--i;
            --j;
        }
        else if (A[i][j] == A[i-1][j] + alpha_gap)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = '-' + b_aligned;
			--i;
        }
        else
        {
            a_aligned = '-' + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned = a[i-1] + a_aligned;
        b_aligned = '-' + b_aligned;
        --i;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned = '-' + a_aligned;
        b_aligned = b[j-1] + b_aligned;
        --j;
    }

    return A[n][m];
}

int CigarAlign(const string &rseq, const string &qseq, const string& cigar, int rbeg, string &r_aligned, string &q_aligned )
{
	r_aligned = rseq.substr( 0, rbeg );
	q_aligned = string( rbeg, '-' );
	char digits[16];
	char symb;
	int npos = 0;
	int ndigit = 0;
	int slen = cigar.size();
	int qbeg = 0;
	while ( npos < slen )
	{
		int symb = cigar[ npos ];
		if ( symb >= '0' && symb <= '9' )
		{
			digits[ ndigit++ ] = symb;
		}
		else
		{
			digits[ ndigit ] = 0;
			ndigit = 0;
			int val = atoi( digits );
			if ( symb == 'D' || symb == 'N' )
			{
				q_aligned += string( val, '-' );
				r_aligned += rseq.substr( rbeg, val );
				rbeg += val;
			}
			else if ( symb == 'I' )
			{
				q_aligned += qseq.substr( qbeg, val );
				r_aligned += string( val, '-' );
				qbeg += val;
			}
			else if ( symb == 'M' || symb == '=' || symb == 'X' )
			{
				q_aligned += qseq.substr( qbeg, val );
				r_aligned += rseq.substr( rbeg, val );
				qbeg += val;
				rbeg += val;
			}
		}
		npos++;
	}
	q_aligned += string( int( rseq.size() ) - rbeg, '-' );
	r_aligned += rseq.substr( rbeg );
}

bool LoadSeq( FILE *sfile, long pos, string& seq )
{
	fseek( sfile, pos, SEEK_SET );
	int state = 0;
	seq.clear();
	do
	{
		int symb = fgetc( sfile );
		if ( state == 0 )
		{
			if ( symb != '>' )
			{
				fprintf( stderr, "assert: '>' expected\n" );
				return false;
			}
			state = 1;
		}
		else if ( state == 1 )
		{
			if ( symb == EOF )
			{
				fprintf( stderr, "sequence expected at the end of file (pos=%ld ftell=%ld)\n", pos, ftell( sfile ) );
				return false;
			}
			if ( symb == '\n' )
			{
				state = 2;
			}
		}
		else if ( state == 2 )
		{
			if ( symb == EOF || symb == '>' ) break;
			if ( symb == '.' ) symb = '-';
			if ( isalpha( symb ) || symb == '-' ) seq.append( 1, toupper( symb ) );
		}
	}
	while ( 1 );
	return true;
}

bool GetHier( const string& id, vector<string>& nkey, bool eukflag )
{
	if ( id.find( ";" ) == -1 )
	{
		nkey.push_back( id );
		return true;
	}
	const int nlevel = 8;
	char tbuf[4096];
	strcpy( tbuf, id.data() );
	vector<string> key;
	key.resize( nlevel );
	const char *kw = "Bacteria;";
	const char *kw2 = "Archaea;";
	const char *kw3 = "Eukaryota;";
	const char *pkw = kw;
	char *pk = strstr( tbuf, kw );
	if ( !pk ) 
	{
		char *pk2 = strstr( tbuf, kw2 );
		if ( !pk2 ) 
		{
			char *pk3 = strstr( tbuf, kw3 );
			if ( !pk3 ) return false;
			else
			{
				pk = pk3;
				pkw = kw3;
				key[0] = kw3;
			}
		}
		else
		{
			pk = pk2;
			pkw = kw2;
			key[0] = "Archaea";
		}
	}
	else key[0] = "Bacteria";
	const char *pk0 = pk;
	pk += strlen( pkw );
	pk += strspn( pk, " \";" );
	char *cpk = pk;
	int lc = 1;
	int len0 = strlen( pk );
	for ( ; lc < nlevel; lc++ )
	{
		int nl = strcspn( cpk, "\";" );
		cpk[ nl ] = 0;
		int clc = lc;
		if ( strncmp( cpk, "otu", 3 ) == 0 ) break;
		key[ clc ] = cpk;
		if ( cpk + strlen( cpk ) + 1 - pk >= len0 ) break;
		char *npl = cpk + strlen( cpk ) + 1;
		npl += strspn( npl, " \";" );
		cpk = npl;
	}
	//nkey.resize( 5 );
	if ( eukflag ) 
	{
		nkey = key;
		int np = strcspn( id.data(), " " ) + 1;
		nkey.back() = string( id, np, pk0 - tbuf - np );
		return true;
	}
	
	if ( key[6].size() > 0 )
	{
		int kc;
		for (  kc = 3; kc < 6; kc++ ) if ( key[kc].substr( key[kc].size() - 3, 3 ) == "les" ) break;
		if ( kc == 6 ) nkey.insert( nkey.begin(), key.begin(), key.end() );
		else
		{
			nkey.insert( nkey.begin(), key.begin(), key.begin() + 2 );
			for ( int kcc = 3; kcc < kc; kcc++ ) nkey.back().append( ";" + key[kcc] );
			nkey.insert( nkey.end(), key.begin() + kc, key.end() );
		}
		while ( nkey.size() > 6 && nkey[6].size() )
		{
			nkey[4].append( ";" + nkey[5] );
			nkey.erase( nkey.begin() + 5 );
		}
		if ( nkey.size() > 6 ) nkey.erase( nkey.begin() + 6, nkey.end() );
	}
	else
	{
		nkey.insert( nkey.begin(), key.begin(), key.begin() + 6 );
	}
	return true;
}

void MergePairs( vector< vector<short> >& ralign, multimap<IPair,IPair>& qpairs, vector<string>& qseqs, vector<string>& qalignseqs )
{
	vector< vector<short> > qalign;
	qalign.resize( ralign.size() );
	int nqseqs = qseqs.size();
	int nrseqs = ralign[0].size();
	for ( int hc = 0; hc < ralign.size(); hc++ )
	{
		qalign[hc].resize( nqseqs, -1 );
		set<IPair> begs;
		for ( int bc = 0; bc < ralign[hc].size(); bc++ )
		{
			if ( ralign[hc][bc] != -1 ) begs.insert( IPair( bc, ralign[hc][bc] ) );
		}
		for ( set<IPair>::iterator it = begs.begin(); it != begs.end(); it++ )
		{
			for ( multimap<IPair,IPair>::iterator it1 = qpairs.lower_bound( *it ); it1 != qpairs.upper_bound( *it ); it1++ )
			{
				int qseqnum = it1->second.first - nrseqs;
				qalign[hc][ qseqnum ] = it1->second.second;
			}
		}
	}
	vector<int> begs( nqseqs, -1 );
	qalignseqs.resize( nqseqs );
	int bound = 0;
	for ( int hc = 0; hc < qalign.size(); hc++ )
	{
		int gmax = 0;
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			if ( qalign[hc][sc] != -1 )
			{
				gmax = max( gmax, bound + qalign[hc][sc] - begs[sc] );
			}
		}
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			if ( qalign[hc][sc] != -1 )
			{
				string sfragm( qseqs[sc], begs[sc] + 1, qalign[hc][sc] - begs[sc] );
				for ( int pc = 0; pc < int( sfragm.size() ) - 1; pc++ ) sfragm[pc] = tolower( sfragm[pc] );
				int gsize = ( gmax - qalignseqs[sc].size()  - sfragm.size() );
				sfragm.insert( 0, gsize, '-' );
				qalignseqs[sc] += sfragm;
				begs[sc] = qalign[hc][sc];
				bound = qalignseqs[sc].size();
			}
		}
	}
	{
		int gmax = 0;
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
				gmax = max( gmax, bound + int( qseqs[sc].size() ) - begs[sc] - 1 );
		}
		for ( int sc = 0; sc < nqseqs; sc++ )
		{
			string sfragm( qseqs[sc], begs[sc] + 1 );
			for ( int pc = 1; pc < int( sfragm.size() ); pc++ ) sfragm[pc] = tolower( sfragm[pc] );
			int gsize = ( gmax - qalignseqs[sc].size()  - sfragm.size() );
			sfragm.append( gsize, '-' );
			qalignseqs[sc] += sfragm;
		}
	}
}

struct RefData
{
	string ref_filename;
	vector<string> ids;
	vector<long> positions;
	map<string,int> index;
	map<int,string> otu;
	map<string,set<int> > rev_otu;
	bool LoadRefs( const char *filename, const char *kwstr );
};

bool RefData::LoadRefs( const char *filename, const char *kwstr )
{
	FILE *ifile = fopen( filename, "rt" );
	if ( !ifile )
	{
		fprintf( stderr, "can't open %s\n", filename );
		return false;
	}
	ref_filename = filename;
	char buf[65535];
	long fpos = 0;
	while ( fgets( buf, 65535, ifile ) )
	{
		strtok( buf, "\r\n" );
		if ( buf[0] == '>' )
		{
			if ( kwstr == 0 || strstr( buf, kwstr ) != 0 )
			{
				index[ string( buf + 1, 0, strcspn( buf + 1, " " ) ) ] = ids.size();
				char *potu = strstr( buf, "otu_" );
				string sotu;
				if ( potu )
				{
					sotu = string( potu, 0, strcspn( potu, " ;,.\n" ) );
				}
				else
				{
					if ( isdigit( buf[1] ) )
					{
						sotu = "otu_" + string( buf + 1, 0, strspn( buf + 1, "0123456789" ) );
					}
				}
				if ( sotu.size() )
				{
					otu[ ids.size() ] = sotu;
					rev_otu[ sotu ].insert( ids.size() );
				}
				ids.push_back( buf + 1 );
				positions.push_back( fpos );
			}
		}
		fpos = ftell( ifile );
	}
	fclose( ifile );
	fprintf( stderr, "input file loaded, index %d, otu %d, rotu %d\n", int( index.size() ), int( otu.size() ), int( rev_otu.size() ) );
	return true;
}

struct SamAlign
{
	string ref;
	int beg;
	string cigar;
};

struct SamData : RefData
{
	map<string,SamAlign> smap;
	map<string, map<int, set<string> > > otusets;
	set<string> rcreads;
	map<string,IPair> otubounds;
	bool LoadSam( const char *sname );
	bool LoadSamGen( const char *lname );
};

IPair parse_cigar( int ref_pos, const char *cigar, string& cleancigar )
{
	char digits[16];
	char symb;
	int npos = 0;
	int ndigit = 0;
	int slen = strlen( cigar );
	bool first = true;
	int beg = 0;
	int length = 0;
	int cbpos = 0;
	int cepos = 0;
	while ( npos < slen )
	{
		int symb = cigar[ npos ];
		if ( symb >= '0' && symb <= '9' )
		{
			digits[ ndigit++ ] = symb;
		}
		else
		{
			digits[ ndigit ] = 0;
			ndigit = 0;
			int val = atoi( digits );
			if ( symb == 'D' || symb == 'N' )
			{
				if ( first ) 
				{
					beg = val;
					cbpos = npos + 1;
				}
				else if ( npos + 1 == slen ) break;
				else 
				{
					length += val;
					cepos = npos;
				}
			}
			else if ( symb == 'M' || symb == '=' || symb == 'X' )
			{
				length += val;
				cepos = npos;
			}
			first = false;
		}
		npos++;
	}
	cleancigar = string( cigar + cbpos, 0, cepos - cbpos + 1 );
	return IPair( beg + ref_pos - 1, beg + length + ref_pos - 1 );
	
}

bool SamData::LoadSam( const char *sname )
{
	FILE *sfile = fopen( sname, "rt" );
	if ( !sfile ) 
	{
		fprintf( stderr, "can't open %s\n", sname );
		return false;
	}
	char buf[65535];
	while ( fgets( buf, 65535, sfile ) )
	{
		if ( buf[0] == '@' ) continue;
		vector<string> fields;
		for ( char *p = strtok( buf, "\t" ); p; p = strtok( 0, "\t" ) ) fields.push_back( p );
		string read_id( fields[0], 0, strcspn( fields[0].data(), " " ) );
		string ref_id( fields[2], 0, strcspn( fields[2].data(), " " ) );
		int flag = atoi( fields[1].data() );
		int ref_pos = atoi( fields[3].data() );
		string& cigar = fields[5];
		string cleancigar;
		IPair cbounds = parse_cigar( ref_pos, cigar.data(), cleancigar );
		if ( index.find( ref_id ) != index.end() ) 
		{
			int ref_ind = index[ ref_id ];
			string otuname = otu[ ref_ind ];
			SamAlign ca;
			ca.ref = ref_id;
			ca.beg = cbounds.first;
			ca.cigar = cleancigar;
			smap[ read_id ] = ca;
			otusets[ otuname ][ ref_ind ].insert( read_id );
			if ( flag & 16 ) rcreads.insert( read_id );
			if ( otubounds.find( otuname ) == otubounds.end() )
			{
				otubounds[ otuname ] = cbounds;
			}
			else
			{
				otubounds[ otuname ] = IPair( min( otubounds[ otuname ].first, cbounds.first ), max( otubounds[ otuname ].second, cbounds.second ) );
			}
		}
	}
	fclose( sfile );
	return true;
}

bool SamData::LoadSamGen( const char *lname )
{
	FILE *lfile = fopen( lname, "rt" );
	if ( !lfile ) 
	{
		fprintf( stderr, "can't open %s\n", lname );
		return false;
	}
	char buf[4096];
	bool first = true;
	while ( fgets( buf, 4096, lfile ) )
	{
		strtok( buf, "\r\n" );
		if ( strlen( buf ) == 0 || buf[0] == '\n' || buf[0] == '\r' ) continue;
		if ( first )
		{
			FILE *tfile = fopen( buf, "rt" );
			if ( !tfile )
			{
				LoadSam( lname );
				break;
			}
			fclose( tfile );
			first = false;
		}
		LoadSam( buf );
	}
	fclose( lfile );
	fprintf( stderr, "sam import done, %d motifs, %d otusets, %d rcreads\n", int( smap.size() ), int( otusets.size() ), int( rcreads.size() ) );
	fflush( stderr );
	return true;
}

struct ReadsData : SamData
{
	map<string,long> read_positions;
	map<string,string> read_files;
	bool LoadReads( const char *fname );
	bool LoadReadsGen( const char *fname );
	set<string> volumes;
};

string samplename( const string& fname )
{
	int pslash = fname.rfind( "/" );
	int pdot = fname.rfind( "." );
	if ( pslash != -1 )
	{
		if ( pdot != -1 && pdot > pslash ) return fname.substr( pslash + 1, pdot - pslash - 1 );
		return fname.substr( pslash + 1 );
	}
	else
	{
		if ( pdot != -1 ) return fname.substr( 0, pdot );
		return fname;
	}
}

bool ReadsData::LoadReads( const char *fname )
{
	FILE *ifile = fopen( fname, "rt" );
	if ( !ifile ) 
	{
		fprintf( stderr, "can't open %s\n", fname );
		return false;
	}
	volumes.insert( samplename( fname ) );
	char buf[65535];
	long fpos = 0;
	while ( fgets( buf, 65535, ifile ) )
	{
		strtok( buf, "\r\n" );
		if ( buf[0] == '>' )
		{
			string id( buf + 1, 0, strcspn( buf + 1, " \t" ) );
			if ( smap.find( id ) != smap.end() )
			{
				read_positions[ id ] = fpos;
				read_files[ id ] = fname;
			}
		}
		fpos = ftell( ifile );
	}
	fclose( ifile );
	return true;
}

bool ReadsData::LoadReadsGen( const char *lname )
{
	FILE *lfile = fopen( lname, "rt" );
	if ( !lfile ) 
	{
		fprintf( stderr, "can't open %s\n", lname );
		return false;
	}
	char buf[4096];
	bool first = true;
	while ( fgets( buf, 4096, lfile ) )
	{
		strtok( buf, "\r\n" );
		if ( strlen( buf ) == 0 || buf[0] == '\n' || buf[0] == '\r' ) continue;
		if ( first )
		{
			FILE *tfile = fopen( buf, "rt" );
			if ( !tfile )
			{
				LoadReads( lname );
				break;
			}
			fclose( tfile );
			first = false;
		}
		LoadReads( buf );
		fprintf( stderr, "file %s, %d /%d reads\n", buf, int( read_positions.size() ), int( read_files.size() )  );
	}
	fclose( lfile );
	fprintf( stderr, "reads import done, %d /%d reads\n", int( read_positions.size() ), int( read_files.size() )  );
	fflush( stderr );
	return true;
}
	
struct MainClass : ReadsData
{
	bool eukflag;
	map<int,int> cref_index;
	vector<int> cref_center;
	vector< map<int,string> > cref_seqs;
	vector< vector< vector<short> > > cref_aligns;
	FILE *of_reads;
	int clcnt1;
	int clcnt2;
	map< vector<string>, map<string,int> > emap;
	
	MainClass() { clcnt1 = clcnt2 = 1; eukflag = false; }
	bool GroupRefs( const string& cotu, double threshold = 0.06, int method = clComplete );
	bool SingleRefs( const string& cotu );
	bool ProcessReads( const string& cotu, double threshold = 0.03, int method = clComplete, int maxgsize = 500 );
};

bool MainClass::GroupRefs( const string& cotu, double threshold, int method )
{
	cref_index.clear();
	cref_center.clear();
	cref_seqs.clear();
	cref_aligns.clear();
	FILE *rfile = fopen( ref_filename.data(), "rt" );
	if ( !rfile )
	{
		fprintf( stderr, "assertion: rfile == 0 at group refs\n" );
		return false;
	}
	vector<string> align;
	vector<int> cindex;
	for ( map<int,set<string> >::iterator it = otusets[ cotu ].begin(); it != otusets[ cotu ].end(); it++ )
	{
		align.resize( align.size() + 1 );
		if ( !LoadSeq( rfile, positions[ it->first ], align.back() ) ) return false;
		cindex.push_back( it->first );
	}
	fclose( rfile );
	vector< vector<float> > distances;
	if ( !CalcDistances( align, distances ) ) return false;
	vector< set<int> > cl_index;
	hierarchicalClustering( distances, threshold, cl_index, method );
	cref_center.resize( cl_index.size() );
	vector<int> ccenter( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		int nmax = 0;
		for ( set<int>::iterator it = cl_index[cc].begin(); it != cl_index[cc].end(); it++ )
		{
			int snum = cindex[ *it ];
			cref_index[ snum ] = cc;
			int ncur = otusets[ cotu ][ snum ].size();
			if ( ncur > nmax )
			{
				cref_center[ cc ] = snum;
				nmax = ncur;
				ccenter[ cc ] = *it;
			}
		}
	}
	int bounds = 20;
	int rbeg = otubounds[ cotu ].first;
	int rend = otubounds[ cotu ].second;
	cref_seqs.resize( cl_index.size() );
	cref_aligns.resize( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		int abeg = -1;
		int aend = -1;
		for ( set<int>::iterator it = cl_index[cc].begin(); it != cl_index[cc].end(); it++ )
		{
			const string& gseq = align[ *it ];
			string refseq;
			for ( int lc = 0; lc < gseq.size(); lc++ )
			{
				if ( refseq.size() == rbeg && ( abeg == -1 || abeg > lc ) )
				{
					abeg = lc;
				}
				if ( refseq.size() == rend && ( aend == -1 || aend < lc ) )
				{
					aend = lc;
				}
				if ( gseq[ lc ] != '-' )
				{
					refseq.append( 1, gseq[ lc ] );
				}
			}
		}
		map<int,string>& gseqs = cref_seqs[cc];
		vector<short> spos( cl_index[cc].size(), -1 );
		for ( int lc = abeg; lc <= aend; lc++ )
		{
			bool fgap = true;
			int c2c = 0;
			for ( set<int>::iterator it = cl_index[cc].begin(); it != cl_index[cc].end(); it++, c2c++ )
			{
				const string& gseq = align[ *it ];
				if ( gseq[ lc ] != '-' )
				{
					spos[c2c] = gseqs[ *it ].size();
					gseqs[ *it ].append( 1, gseq[lc] );
					fgap = false;
				}
			}
			if ( fgap ) continue;
			cref_aligns[cc].push_back( spos );
		}
	}
	return true;
}

bool MainClass::SingleRefs( const string& cotu )
{
	cref_index.clear();
	cref_center.clear();
	cref_seqs.clear();
	cref_aligns.clear();
	FILE *rfile = fopen( ref_filename.data(), "rt" );
	if ( !rfile )
	{
		fprintf( stderr, "assertion: rfile == 0 at group refs\n" );
		return false;
	}
	vector<string> rseqs;
	vector< set<int> > cl_index;
	for ( map<int,set<string> >::iterator it = otusets[ cotu ].begin(); it != otusets[ cotu ].end(); it++ )
	{
		rseqs.resize( rseqs.size() + 1 );
		if ( !LoadSeq( rfile, positions[ it->first ], rseqs.back() ) ) return false;
		set<int> singleset;
		singleset.insert( it->first );
		cl_index.push_back( singleset );
	}
	fclose( rfile );
	cref_center.resize( cl_index.size() );
	vector<int> ccenter( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		cref_center[ cc ] = *( cl_index[cc].begin() );
		cref_index[ *( cl_index[cc].begin() ) ] = cc;
	}
	int bounds = 20;
	int rbeg = otubounds[ cotu ].first;
	int rend = otubounds[ cotu ].second;
	cref_seqs.resize( cl_index.size() );
	cref_aligns.resize( cl_index.size() );
	for ( int cc = 0; cc < cl_index.size(); cc++ )
	{
		int abeg = rbeg;
		int aend = rend;
		for ( int lc = abeg; lc <= aend; lc++ )
		{
			vector<short> spos( 1, lc - abeg );
			cref_aligns[cc].push_back( spos );
		}
		cref_seqs[cc][cc] = rseqs[cc].substr( abeg, aend - abeg + 2 );
	}
	return true;
}

string itostr( int n )
{
	char buf[16];
	sprintf( buf, "%d", n );
	return string( buf );
}

bool MainClass::ProcessReads( const string& cotu, double threshold, int clmethod, int maxgsize )
{
    int alpha[alphabets][alphabets];
    for (size_t i = 0; i < alphabets; ++i)
    {
        for (size_t j = 0; j < alphabets; ++j)
        {
            if (i == j) alpha[i][j] = 0;
            else alpha[i][j] = 1;
        }
    }
    int gap_penalty = 2;

	FILE *rfile = 0;
	string rfname;
	string sfname;
	int seqc = 0;
	for ( map< int,set<string> >::iterator it0 = otusets[ cotu ].begin(); it0 != otusets[ cotu ].end(); it0++, seqc++ )
	{
		fprintf( stderr, "  ref %d, %d reads\n", seqc, int( it0->second.size() ) );
		vector< vector<string> > rseqs;
		vector< multimap<IPair,IPair> > ralign;
		vector< vector<double> > rscores;
		vector< vector<string> > rsamples;
		vector< vector<string> > rids;
		rseqs.resize( cref_seqs.size() );
		ralign.resize( cref_seqs.size() );
		rscores.resize( cref_seqs.size() );
		rsamples.resize( cref_seqs.size() );
		rids.resize( cref_seqs.size() );
		for ( set<string>::iterator it = it0->second.begin(); it != it0->second.end(); it++ )
		{
			if ( read_files.find( *it ) == read_files.end() ) continue;
			if ( rfile == 0 || rfname != read_files[ *it ] )
			{
				if ( rfile ) fclose( rfile );
				rfname = read_files[ *it ];
				rfile = fopen( rfname.data(), "rt" );
				if ( !rfile )
				{
					fprintf( stderr, "can't open %s [%s]\n", rfname.data(), it->data() );
					return false;
				}
				sfname = samplename( rfname );
			}
			string rseq;
			if ( !LoadSeq( rfile, read_positions[ *it ], rseq ) ) return false;
			if ( rcreads.find( *it ) != rcreads.end() )
			{
				string rcseq;
				for ( int lc = 0; lc < rseq.size(); lc++ )
				{
					int clet = rseq[ rseq.size() - lc - 1 ];
					int nlet;
					switch ( clet )
					{
						case 'A': nlet = 'T'; break;
						case 'T': nlet = 'A'; break;
						case 'G': nlet = 'C'; break;
						case 'C': nlet = 'G'; break;
						nlet = 'N';
					}
					rcseq.append( 1, nlet );
				}
				rseq = rcseq;
			}
			int cluster = cref_index[ it0->first ];
			string qalign;
			string salign;
			string& sseq = cref_seqs[cluster][seqc];
			//int al_res = NWAlign( rseq, sseq, gap_penalty, alpha, qalign, salign );
			CigarAlign( sseq, rseq, smap[ *it ].cigar, smap[ *it ].beg - otubounds[ cotu ].first, salign, qalign );
			int nseq = rseqs[ cluster ].size();
			rseqs[ cluster ].push_back( rseq );
			rsamples[ cluster ].push_back( sfname );
			rids[ cluster ].push_back( *it );
			int b1 = 0;
			int b2 = 0;
			int nident = 0;
			int cref_size = cref_aligns[ cluster ][0].size(); 
			for ( int lc = 0; lc < qalign.size(); lc++ )
			{
				if ( qalign[lc] == '-' ) b2++;
				else if ( salign[lc] == '-' ) b1++;
				else
				{
					IPair p1( seqc, b2 );
					IPair p2( nseq + cref_size, b1 );
					ralign[ cluster ].insert( pair<IPair,IPair>( p1, p2 ) );
					ralign[ cluster ].insert( pair<IPair,IPair>( p2, p1 ) );
					if ( qalign[ lc ] == salign[ lc ] ) nident++;
					b1++;
					b2++;
				}
			}
			//fprintf( stderr, "%s\n%s\nident %d qseq %d sseq %d al_res %d read %s bounds %d %d\n", qalign.data(), salign.data(), nident, int( rseq.size() ), int( sseq.size() ), rcreads.find( *it ) != rcreads.end() ? 1 : 0, it->data(), otubounds[ cotu ].first, otubounds[ cotu ].second );
			rscores[ cluster ].push_back( double( nident ) / rseq.size() ); 
		}
		for ( int cc = 0; cc < cref_seqs.size(); cc++ )
		{
			int otucl = clcnt1;
			clcnt1++;
			vector<string> hkey;
			GetHier( ids[ cref_center[cc] ], hkey, eukflag );
			hkey.push_back( itostr( otucl ) );
			if ( rseqs[cc].size() > maxgsize )
			{
				int rcl = clcnt2;
				clcnt2++;
				double sumscores = 0;
				for ( int sc = 0; sc < rseqs[cc].size(); sc++ )
				{
					sumscores += rscores[cc][sc];
				}
				double pident = 100. * sumscores / rseqs[cc].size();
				char nbuf[16];
				sprintf( nbuf, "%d_%4.1f", rcl, pident );
				vector<string> emkey = hkey;
				emkey.push_back( nbuf );
				for ( int sc = 0; sc < rseqs[cc].size(); sc++ )
				{
					fprintf( of_reads, ">%s %s %s %d_%d_%4.1f %4.1f\n%s\n", rids[cc][sc].data(), rsamples[cc][sc].data(), cotu.data(), otucl, rcl, pident, 100. * rscores[cc][sc], rseqs[cc][sc].data() );
					emap[ emkey ][ rsamples[cc][sc] ]++;
				}
				fflush( of_reads );
				continue;
			}
			fprintf( stderr, "   %d reads to cluster\n", int( rseqs[cc].size() ) );
			vector<string> ralignseqs;
			MergePairs( cref_aligns[cc], ralign[cc], rseqs[cc], ralignseqs );
			vector< vector<float> > distances;
			if ( !CalcDistances( ralignseqs, distances ) ) return false;
			vector< set<int> > cl_index;
			hierarchicalClustering( distances, threshold, cl_index, clmethod );
			for ( int c2c = 0; c2c < cl_index.size(); c2c++ )
			{
				int rcl = clcnt2;
				clcnt2++;
				double sumscores = 0;
				for ( set<int>::iterator it = cl_index[c2c].begin(); it != cl_index[c2c].end(); it++ )
				{
					sumscores += rscores[cc][ *it ];
				}
				double pident = 100. * sumscores / cl_index[c2c].size();
				char nbuf[16];
				sprintf( nbuf, "%d_%4.1f", rcl, pident );
				vector<string> emkey = hkey;
				emkey.push_back( nbuf );
				for ( set<int>::iterator it = cl_index[c2c].begin(); it != cl_index[c2c].end(); it++ )
				{
					int sc = *it;
//					fprintf( of_reads, ">%s %s %s %d_%d_%4.1f %4.1f\n%s\n", rids[cc][sc].data(), rsamples[cc][sc].data(), cotu.data(), otucl, rcl, pident, 100. * rscores[cc][sc], rseqs[cc][sc].data() );
					fprintf( of_reads, ">%s %s %s %d_%d_%4.1f %4.1f\n%s\n", rids[cc][sc].data(), rsamples[cc][sc].data(), cotu.data(), otucl, rcl, pident, 100. * rscores[cc][sc], ralignseqs[sc].data() );
					emap[ emkey ][ rsamples[cc][sc] ]++;
				}
			}
			fflush( of_reads );
		}
	}
	if ( rfile ) fclose( rfile );
	return true;
}

int main( int argc, char **argv )
{
	if ( argc == 1 )
	{
		fprintf( stderr, "arguments: ref_fasta sam/list reads/list out.fa out_table.txt\n" );
		return 1;
	}
	const char *refname = argv[1];
	const char *samname = argv[2];
	const char *readsname = argv[3];
	const char *of1_name = argv[4];
	const char *of2_name = argv[5];
	const char *kwstr = 0;
	double thresh1 = 0.6;
	double thresh2 = 0.03;
	int clmethod1 = clComplete;
	int clmethod2 = clComplete;
	int maxgsize = 500;
	
	MainClass data;
	if ( !data.LoadRefs( refname, kwstr ) ) return 1;
	if ( !data.LoadSamGen( samname ) ) return 1;
	if ( !data.LoadReadsGen( readsname ) ) return 1;
	data.of_reads = fopen( of1_name, "wt" );
	for ( map<string, map<int,set<string> > >::iterator it = data.otusets.begin(); it != data.otusets.end(); it++ )
	{
		fprintf( stderr, "%s: %d refs... ", it->first.data(), int( it->second.size() ) ); 
		fflush( stderr );
		//if ( !data.GroupRefs( it->first, thresh1, clmethod1 ) ) return 1;
		if ( !data.SingleRefs( it->first ) ) return 1;
		fprintf( stderr, "%d ref clusters\n", int( data.cref_seqs.size() ) );
		if ( !data.ProcessReads( it->first, thresh2, clmethod2, maxgsize ) ) return 1;
	}
	int nlevels = data.emap.begin()->first.size();
	FILE *of_table = fopen( of2_name, "wt" );
	for ( int lc = 1; lc <= nlevels; lc++ ) fprintf( of_table, "%d\t", lc );
	for ( set<string>::iterator it = data.volumes.begin(); it != data.volumes.end(); it++ )
	{
		fprintf( of_table, "%s\t", it->data() );
	}
	fprintf( of_table, "\n" );
	for ( map< vector<string>, map<string,int> >::iterator it = data.emap.begin(); it != data.emap.end(); it++ )
	{
		for ( int lc = 0; lc < nlevels; lc++ ) fprintf( of_table, "%s\t", it->first[lc].data() );
		for ( set<string>::iterator it1 = data.volumes.begin(); it1 != data.volumes.end(); it1++ )
		{
			int count = 0;
			if ( it->second.find( *it1 ) != it->second.end() ) count = it->second[ *it1];
			fprintf( of_table, "%d\t", count );
		}
		fprintf( of_table, "\n" );
	}
	fclose( of_table );
	return 0;
}
	
	

	
	