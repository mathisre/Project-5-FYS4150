#ifndef IO_H
#define IO_H
#include <fstream>
class System;
using std::ofstream;

class IO
{
private:
    ofstream file;
public:
    IO(std::string filename);
    ~IO();

    void saveState(System &system);
    void open(std::string filename);
    void close();

};
#endif
