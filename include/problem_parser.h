#ifndef PROBLEM_FILE_PARSER
#define PROBLEM_FILE_PARSER

#include <fstream>
#include <string>
#include <utility>
#include "Matrix.h"


using std::string;
using namespace my_matrix;

/*
 * Класс основан TSPLIB 95 инфа о библе http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
 * Не обращайте внимание на заглавные имена полей класса, сделанно для удобство парсинга
 * При парсинге не забывайте закрывать файл
 * */

class problem_parser {
protected:
    /* Внутренние поля */
    bool started_read = false;
    std::ifstream file_point_; /* Указатель на каретку в файле */
    string curent_line_; /* Текущая строка */
public:
    /* Поля файла */
    string filename_; /* Имя файла */
    string problem_name_;
    string task_type_; /* Тип задания */
    string comment_;
    size_t dimension_{}; /* Размерность задачи */
    size_t capacity_{};  /* Указывает вместимость грузовика в CVRP (задача с грузовиком)*/
    string edge_weight_type_; /*Тип весов узлов 2d евклид, или явная задача весов*/
    string edge_weight_format_; /* Формат весов если заданы явно */
    string edge_data_format_; /* Описывает формат вершин графа, если граф не полный*/
    string node_coord_type_; /* Указывает, связаны ли координаты с каждым узлом */
    string display_data_type_; /* Формат весов */
    /* Секция данных */
    int data_sector_type_ = 0;
protected:
    std::ifstream StartParse(const string& filename = "");
    static int CheckIfDataSector(const string& str);
    int FindNextDataSector();
    static double SplitLine(std::stringstream& ss);
public:
    problem_parser();
    explicit problem_parser(string  filename): filename_(std::move(filename)){
        file_point_ = StartParse();
        started_read = true;
    };
    Matrix* MatrixParseFile(const string& filename = "");
};

#endif