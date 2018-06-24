
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <string>
#include <vector>
#include <map>
#include <set>

using namespace std;

struct Record
{
	string aseq;
	string seq;
	string id;
	string species;
};

const int maxphsize = 2000;

void SelectPhyla( const char *iname, map<string,int>& phyla )
{
	FILE *ifile = fopen( iname, "rt" );
	char buf[100000];
	string id;
	string tax;
	bool first = true;
	while ( fgets( buf, 100000, ifile ) )
	{
		if ( strncmp( buf, "\tspecies", 8 ) == 0 )
		{
			if ( first )
			{
				first = false;
			}
			else if ( tax.find( "Bacteria" ) != -1 && id.size() )
			{
				vector<string> etax;
				char tbuf[4096];
				strcpy( tbuf, tax.data() );
				char *p1 = strchr( tbuf, ';' );
				if ( p1 )
				{
					char *p2 = strchr( p1 + 1, ';' );
					if ( p2 )
					{
						char *p3 = strchr( p2 + 1, ';' );
						if ( p3 )
						{
							*p3 = 0;
						}
					}
					phyla[ tbuf ]++;
				}
			}
			id.clear();
			tax.clear();
		}
		if ( strncmp( buf, "\t\tacc", 5 ) == 0 )
		{
			char *p = strtok( buf + 5, "\t\"\n" );
			if ( p ) id = p;
		}
		if ( strncmp( buf, "\t\ttax_slv", 9 ) == 0 )
		{
			char *p = strtok( buf + 9, "\t\"\n" );
			if ( p ) tax = p;
		}
	}
	fclose( ifile );
}

void SelectPhylaKW( const char *iname, set<string>& keywords, set<string>& ukeywords, map<string,int>& phyla )
{
	FILE *ifile = fopen( iname, "rt" );
	char buf[100000];
	string id;
	string tax;
	bool first = true;
	while ( fgets( buf, 100000, ifile ) )
	{
		if ( strncmp( buf, "\tspecies", 8 ) == 0 )
		{
			if ( first )
			{
				first = false;
			}
			else if ( tax.find( "Bacteria" ) != -1 && id.size() )
			{
				string ckw;
				for ( set<string>::iterator it = keywords.begin(); it != keywords.end(); it++ )
				{
					if ( strncmp( it->data(), tax.data(), it->size() ) == 0 )
					{
						ckw = *it;
						break;
					}
				}
				if ( ckw.size() > 0 )
				{
					vector<string> etax;
					char tbuf[4096];
					strcpy( tbuf, tax.data() );
					char *p1 = strchr( tbuf + ckw.size() + 1, ';' );
					if ( p1 )
					{
						char *p2 = strchr( p1 + 1, ';' );
						if ( p2 )
						{
							*p2 = 0;
						}
						phyla[ tbuf ]++;
						ukeywords.insert( ckw );
					}
				}
			}
			id.clear();
			tax.clear();
		}
		if ( strncmp( buf, "\t\tacc", 5 ) == 0 )
		{
			char *p = strtok( buf + 5, "\t\"\n" );
			if ( p ) id = p;
		}
		if ( strncmp( buf, "\t\ttax_slv", 9 ) == 0 )
		{
			char *p = strtok( buf + 9, "\t\"\n" );
			if ( p ) tax = p;
		}
	}
	fclose( ifile );
}


void ProcessPhyla( const char *iname, const string& phyla, long *tucnt, FILE *of1, FILE *of2 )
{
	FILE *ifile = fopen( iname, "rt" );
	char buf[100000];
	string id;
	string tax;
	string tax2;
	string species;
	string seq;
	bool first = true;
	long rc = 0;
	map<vector<string>,vector<Record> > gdata;
	
	while ( fgets( buf, 100000, ifile ) )
	{
		if ( strncmp( buf, "\tspecies", 8 ) == 0 )
		{
			if ( first )
			{
				first = false;
			}
			else if ( tax.find( "Bacteria" ) != -1 && id.size() && strncmp( tax.data(), phyla.data(), phyla.size() ) == 0 && rc < maxphsize )
			{
				vector<string> etax;
				char tbuf[4096];
				strcpy( tbuf, tax.data() );
				for ( char *p = strtok( tbuf, ";" ); p; p = strtok( 0, ";" ) ) etax.push_back( p );
				Record r;
				r.aseq = seq;
				for ( int lc = 0; lc < seq.size(); lc++ ) if ( seq[lc] != '-' && seq[lc] != '.' ) r.seq.append( 1, seq[lc] );
				r.id = id;
				if ( tax2.find( "Bacteria" ) != -1 && species.size() && ( species.find( "uncultured" ) == -1 ) )
				{
					r.species = species;
				}
				vector<string> key;
				key.insert( key.begin(), etax.begin(), etax.begin() + min( 6, int( etax.size() ) ) );
				gdata[key].push_back( r );
				rc++;
				fprintf( stderr, "%d %s %s\n", rc, tax.data(), id.data() );
				fflush( stderr );
			}
			id.clear();
			tax.clear();
			seq.clear();
			tax2.clear();
			species.clear();
		}
		if ( strncmp( buf, "\t\tacc", 5 ) == 0 )
		{
			char *p = strtok( buf + 5, "\t\"\n" );
			if ( p ) id = p;
		}
		if ( strncmp( buf, "\t\ttax_slv", 9 ) == 0 )
		{
			char *p = strtok( buf + 9, "\t\"\n" );
			if ( p ) tax = p;
		}
		if ( strncmp( buf, "\t\ttax_embl\t", 11 ) == 0 )
		{
			char *p = strtok( buf + 11, "\t\"\n" );
			if ( p ) tax2 = p;
		}
		if ( strncmp( buf, "\t\tfull_name", 11 ) == 0 )
		{
			char *p = strtok( buf + 11, "\t\"\n" );
			if ( p ) species = p;
		}
		if ( strncmp( buf, "\t\t\tdata", 7 ) == 0 )
		{
			char *p = strchr( buf, '"' );
			if ( p )
			{
				for ( ; *p > 0; p++ ) if ( isalpha( *p ) || *p == '-' || *p == '.' ) seq.append( 1, *p );
			}
		}
	}
	fclose( ifile );
	fprintf( stderr, "%s %d entries %ld records\n", phyla.data(), int( gdata.size() ), rc );
	fflush( stderr );
	int ocnt = 1;
	int otucnt = 1;
	const char *tuname[] = { "kitu", "phtu", "cltu", "ortu", "fatu", "getu", "sptu" };
	const char *tulevel[] = { "0", "0.5", "0.7", "0.78", "0.85", "0.92", "0.97" };
	for ( map<vector<string>,vector<Record> >::iterator it = gdata.begin(); it != gdata.end(); it++ )
	{
		string otax;
		for ( int tc = 0; tc < it->first.size(); tc++ ) otax.append( it->first[tc] + "; " );
		vector<Record> *vr = &( it->second );
		vector<Record> nr;
		int levelc = 0;
		int nlevel = 6 - it->first.size();
		vector<int> rcenters;
		vector< vector<string> > rtax;
		vector< set<int> > rsets;
		do
		{
			FILE *ucfile = fopen( "uctmp.fa", "wt" );
			for ( int rc = 0; rc < vr->size(); rc++ )
			{
				fprintf( ucfile, ">%d\n%s\n", rc, ( *vr )[rc].seq.data() );
			}
			fclose( ucfile );
			int r1 = system( "/mnt/disk3/sergey/miniconda2/pkgs/qiime-1.9.1-np110py27_0/bin/uclust --sort uctmp.fa --output uctmp_s.fa" );
			int r2 = system( ( string( "/mnt/disk3/sergey/miniconda2/pkgs/qiime-1.9.1-np110py27_0/bin/uclust --input uctmp_s.fa --uc tmp.uc --id " ) + tulevel[ 6 - levelc ] ).data() );
			map<int,int> clcenters;
			map<int,string> clspecies;
			map<int,set<int> > clsets;
			FILE *urfile = fopen( "tmp.uc", "rt" );
			if ( !urfile ) continue;
			while ( fgets( buf, 4096, urfile ) )
			{
				vector<string> fields;
				for ( char *p = strtok( buf, "\t\n" ); p; p = strtok( 0, "\t\n" ) ) fields.push_back( p );
				if ( fields.size() < 9 ) continue;
				string id = fields[8];//.substr( 0, strcspn( fields[8].data(), " " ) );
				if ( fields[0] == "S" ) 
				{
					int rn = atoi( id.data() );
					int cn = atoi( fields[1].data() );
					if ( rn >= 0 && rn < it->second.size() )
					{
						Record& cr = ( *vr )[rn];
						clcenters[ cn ] = rn;
						if ( cr.species.size() )
						{
							clspecies[ cn ] = cr.species;
						}
						if ( levelc > 0 )
						{
							clsets[ cn ].insert( rsets[ rn ].begin(), rsets[ rn ].end() );
						}
					}
				}
				if ( fields[0] == "H" )
				{
					int rn = atoi( id.data() );
					int cn = atoi( fields[1].data() );
					if ( rn >= 0 && rn < it->second.size() && clcenters.find( cn ) != clcenters.end() )
					{
						Record& cr = ( *vr )[rn];
						if ( cr.species.size() && clspecies.find( cn ) == clspecies.end() )
						{
							clcenters[ cn ] = rn;
							clspecies[ cn ] = cr.species;
						}
						if ( levelc > 0 )
						{
							clsets[ cn ].insert( rsets[ rn ].begin(), rsets[ rn ].end() );
						}
					}
				}
			}
			fclose( urfile );
			if ( levelc == 0 )
			{
				for ( map<int,int>::iterator it1 = clcenters.begin(); it1 != clcenters.end(); it1++ )
				{
					rcenters.push_back( it1->second );
					if ( nlevel > 0 )
					{
						nr.push_back( it->second[ it1->second ] );
						set<int> sset;
						sset.insert( rcenters.size() - 1 );
						rsets.push_back( sset );
					}
				}
				rtax.resize( rcenters.size() );
			}
			else
			{
				rsets.clear();
				vector<Record> pnr = nr;
				nr.clear();
				for ( map<int,int>::iterator it1 = clcenters.begin(); it1 != clcenters.end(); it1++ )
				{
					char nbuf[32];
					sprintf( nbuf, "%s_%ld", tuname[ 6 - levelc ], tucnt[ 6 - levelc ] );
					tucnt[ 6 - levelc ]++;
					for ( set<int>::iterator it2 = clsets[ it1->first ].begin(); it2 != clsets[ it1->first ].end(); it2++ )
					{
						rtax[ *it2 ].push_back( nbuf );
					}
					if ( nlevel > levelc )
					{
						nr.push_back( pnr[ it1->second ] );
						rsets.push_back( clsets[ it1->first ] );
					}
				}
			}
			if ( levelc >= nlevel )
			{
				break;
			}
			levelc++;
			vr = &(nr);
		}
		while ( 1 );
		map<string,int> pgenus;
		for ( int sc = 0; sc < rcenters.size(); sc++ )
		{
			string ctax = otax;
			for ( int tc = int( rtax[sc].size() ) - 1; tc >= 0; tc-- ) ctax.append( rtax[sc][tc] + "; " );
			int otunum = tucnt[ 6 ];
			if ( rtax[sc].size() )
			{
				if ( pgenus.find( rtax[sc][0] ) == pgenus.end() )
				{
					pgenus[ rtax[sc][0] ] = tucnt[ 6 ]++;
				}
				otunum = pgenus[ rtax[sc][0] ];
			}
			Record& cr = it->second[ rcenters[ sc ] ];
			species = "Unk_sp";
			if ( cr.species.size() ) species = cr.species;
			char oid[4096];
			sprintf( oid, "%ld %s %s %s otu_%d", tucnt[7], cr.id.data(), species.data(), ctax.data(), otunum );
			fprintf( of1, ">%s\n%s\n", oid, cr.aseq.data() );
			fprintf( of2, ">%s\n%s\n", oid, cr.seq.data() );
			fflush( of1 );
			fflush( of2 );
			tucnt[7]++;
		}
		tucnt[ 6 ] ++;
	}
}


int main( int argc, char **argv )
{
	if ( argc == 1 )
	{
		fprintf( stderr, "arguments: silva oprefix\n" );
		return 1;
	}
	const char *iname = argv[1];
	map<string,int> nphyla;
	SelectPhyla( iname, nphyla );
	map<string,int> phyla;
	for ( int ic = 0; ic < 3; ic++ )
	{
		set<string> keywords;
		for ( map<string,int>::iterator it = nphyla.begin(); it != nphyla.end(); )
		{
			if ( it->second < maxphsize )
			{
				phyla.insert( *it );
				nphyla.erase( it++ );
			}
			else
			{
				keywords.insert( it->first );
				it++;
			}
		}
		set<string> ukeywords;
		map<string,int> pphyla;
		SelectPhylaKW( iname, keywords, ukeywords, pphyla );
		for ( map<string,int>::iterator it = nphyla.begin(); it != nphyla.end(); it++ )
		{
			if ( ukeywords.find( it->first ) == ukeywords.end() )
			{
				phyla.insert( *it );
			}
		}
		nphyla.clear();
		for ( map<string,int>::iterator it = pphyla.begin(); it != pphyla.end(); it++ )
		{
			if ( it->second < maxphsize )
			{
				phyla.insert( *it );
			}
			else
			{
				nphyla.insert( *it );
			}
		}
		if ( nphyla.size() == 0 ) break;
		fprintf( stderr, "ic %d phyla %d nphyla %d\n", ic, int(  phyla.size() ), int( nphyla.size() ) );
	}
				
	fprintf( stderr, "%d phyla\n", int( phyla.size() ) );
	fflush( stderr );
	FILE *of1 = fopen( ( string( argv[2] ) + "_a.fa" ).data(), "wt" );
	FILE *of2 = fopen( ( string( argv[2] ) + "_c.fa" ).data(), "wt" );
	long tucnt[8] = { 0 };
	for ( map<string,int>::iterator it = phyla.begin(); it != phyla.end(); it++ )
	{
		fprintf( stderr, "Phyla %s %d\n", it->first.data(), it->second );
		fflush( stderr );
		ProcessPhyla( iname, it->first, tucnt, of1, of2 );
	}
	fclose( of1 );
	fclose( of2 );
	return 0;
}
