#define TEST_MODE
#define USE_ARRAY
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include "array.hpp"
/*	some MATLAB code to get data out...
fid = fopen('binfile.bin', 'r');
fseek(fid, -64, 'eof');

A = fread(fid, 64, '*char');



fclose(fid);
*/
using namespace std;


template<class T>
inline string to_string(const T& t) {
	stringstream ss;
	ss << t;
	return ss.str();
}




/************************************************************/

class File {
public:
	string filename;
	struct MetadataFooter{
		char creatorVersion[256];
		char creatorApplication[256];
		char creatorDate[256];
		char creatorSystem[256];
		char creatorComputer[256];
	};
	
	struct FileFooter
	{
		unsigned long int metadataFooterSize;  // = sizeof(MetadataFooter)
		char magicString[256];   // a unique identifier for the format: maybe "MYFILEFMT"
	};
	
	File(){SetMetaData();};
	File(string name) : filename(name) {SetMetaData();};

	void Save() {
		//	Open file to write
		vOpenWrite();
		//	Write content of buffer
		WriteFromBuffer();
		//	Write MetaData
		vWriteMetaData();
		//	Close file
		file.close();
		//	Clear output buffer
		ClearOutputBuffer();
	};

	void Read() {
		//	Open file to read
		vOpenRead();
		//	Clear input buffer
		ClearInputBuffer();
		//	Write content to buffer
		ReadToBuffer();
		//	Close file
		file.close();
	};

protected:
	fstream file;
	stringstream inBuffer, outBuffer;
	void WriteFromBuffer(){file << outBuffer.str();};
	void ReadToBuffer(){inBuffer << file.rdbuf();}
	void ClearOutputBuffer(){outBuffer.str("");};
	void ClearInputBuffer(){inBuffer.str("");};
	virtual void vOpenRead() = 0;
	virtual void vOpenWrite() = 0;
	//virtual void vWriteMetaData() = 0;
	void vWriteMetaData() {
		file << MetaData.creatorVersion << endl;
		file << MetaData.creatorApplication << endl;
		file << MetaData.creatorDate << endl;
		file << MetaData.creatorSystem << endl;
		file << MetaData.creatorComputer << endl;
		file << Footer.metadataFooterSize << endl;
		file << Footer.magicString << endl;		
		};
	void SetMetaData(){
		strcpy (MetaData.creatorVersion,"Version 1 Release 0");
		strcpy (MetaData.creatorApplication,"Combined Vorticity Transport Panel Method Code");
		strcpy (MetaData.creatorDate,"12/34/23");
		strcpy (MetaData.creatorSystem,"Linux x64");
		strcpy (MetaData.creatorComputer,"Chameleon");
		strcpy (Footer.magicString,"CHAMELEON");
		Footer.metadataFooterSize = sizeof(MetadataFooter);
	}

	MetadataFooter MetaData;
	FileFooter Footer;

	enum exception {
		NO_FILE
		};
	};

/************************************************************/
	class binaryFile: public File {
	public:

		binaryFile() : File() {};
		binaryFile(string name) : File(name) {};

			template <typename T> void Write(T &in) {	//	would it be more efficient to put this in as a pointer?		
		outBuffer.write((char*) &in, sizeof(in));
		cout << "Writing" << endl;
	};

			template <typename T> void ReCast(T &in) {};

private:
	void vOpenWrite() {file.open(filename.c_str(), ios::out | ios::binary | ios::app);};
	void vOpenRead() {file.open(filename.c_str(), ios::in | ios::binary);};
	//void vWriteMetaData() {};
};
/************************************************************/
class asciiFile: public File {
public:
	asciiFile() : File() {};
	asciiFile(string name) : File(name) {};

	template <typename T> void Write(T in) {
	outBuffer << in << "\n";
};

string str()
{
	Read();
	return inBuffer.str();
}

private:
	void vOpenWrite() {file.open(filename.c_str(), ios::out | ios::trunc);};
	void vOpenRead() {file.open(filename.c_str());};
};
/************************************************************/
int main () {
	asciiFile asciiF(string("test.txt"));
	binaryFile binaryF(string("binfile.bin"));
	cout << asciiF.filename << " " << binaryF.filename << endl;
	int a = 3;
	double b = 1231214.12412;

	Array <double> c(3);
	c[0] = 123.456;
	c[1] = 234.567;
	c[2] = -9876543.21;

	asciiF.Write(a);
	asciiF.Write(b);
	for (int i = 0; i < 12; ++i)
		asciiF.Write(i);

	asciiF.Write(string("Willies"));
	asciiF.Save();

	asciiF.Read();
	cout << asciiF.str();


	binaryF.Write(a);
	binaryF.Write(a);
	binaryF.Write(b);

	binaryF.Save();

	binaryF.Read();


	return 0;
}