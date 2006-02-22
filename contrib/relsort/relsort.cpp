#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <iostream>
#include <iomanip>
#include <sstream>

#include <glob.h>
#include <map>
using namespace std;

const off_t limit = 1ul << 31;

typedef multimap<off_t,char*> fmap_t;

int main(int argc, char* argv[]) {

    if(argc<3) {
	cout << "Usage: " << argv[0] << " [-d] <ibase> <obase>" << endl;
        return 1;
    }

    bool delsrc = true;
    string ibase(argv[1]); ibase += "*";
    string obase(argv[2]);

    off_t totalsize = 0;

    fmap_t M;

    glob_t glb;
    glob(ibase.c_str(), GLOB_NOSORT, 0, &glb);

//    clog << "Quering size ";
    for(int i=0;i<glb.gl_pathc;++i) {
	//cout << glb.gl_pathv[i] << endl;
//        clog << "."; clog.flush();

        struct stat st;
        stat(glb.gl_pathv[i], &st);

	M.insert(make_pair(st.st_size,glb.gl_pathv[i]));
        totalsize += st.st_size;
    }
    clog << endl << "Files: " << glb.gl_pathc << "\tTotal size: " << totalsize << endl;

    int blk = 0;

    while(!M.empty()) {
	off_t sz = limit;
	clog << "Block " << blk << " : ";
        cout << "cat ";

	for(fmap_t::reverse_iterator im=M.rbegin();im!=M.rend();) {
	    if( im->first <= sz ) {
		sz -= im->first;
		cout << im->second;

		if(++im == M.rend()) {
		    M.erase(im.base());
		    break;
		}
		cout << " ";
      	        M.erase(im.base());
	    }
	    else ++im;
	}
	clog << "size " << limit - sz << endl;
        clog.flush();

	ostringstream of;
        of << " " << obase << "." << setfill('0') << setw(3) << blk;
	cout << ">" << of.str() << endl;
        ++blk;
    }

    globfree(&glb);
    return 0;
}
