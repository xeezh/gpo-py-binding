#include "../../include/problem_parser.h"
#include "MathCoordinate.cpp"
#include <sstream>

problem_parser::problem_parser() = default;
std::ifstream problem_parser::StartParse(const string& filename) {
    if (!filename.empty()) {
        filename_ = filename;
    } else if (filename_.empty()) {
        throw std::invalid_argument("Filename is not provided");
    }
    std::ifstream in_file(filename_);
    if (!in_file) {
        throw std::invalid_argument("Error when file open");
    }

    string delimiter = ":";
    std::string name;
    std::string value;
    size_t pos;
    size_t erase_pos;

    bool is_data_segment = false; /* Указывает расположена ли каретка на сегменте с данными */
    while (!is_data_segment) {
        std::getline(in_file, curent_line_);
        string line = curent_line_;
        erase_pos = line.find('\r');
        while (erase_pos != std::string::npos) {
            line.erase(erase_pos, 1);
            erase_pos = line.find('\r');
        }
        erase_pos = line.find(' ');
        while (erase_pos != std::string::npos) {
            line.erase(erase_pos, 1);
            erase_pos = line.find(' ');
        }


        pos = line.find(delimiter);

        if (pos != std::string::npos) {
            name = line.substr(0, pos);
            value = line.substr(pos + delimiter.length());
        } else{
            name = line;
        }



        if (name == "NAME") {
            problem_name_ = value;
        } else if (name == "TYPE") {
            task_type_ = value;
        } else if (name == "COMMENT") {
            pos = curent_line_.find(delimiter);
            value = curent_line_.substr(pos + delimiter.length());
            comment_ += value;
        } else if (name == "DIMENSION") {
            dimension_ = std::stoi(value);
        } else if (name == "CAPACITY") {
            capacity_ = std::stoi(value);
        } else if (name == "EDGE_WEIGHT_TYPE") {
            edge_weight_type_ = value;
        } else if (name == "EDGE_WEIGHT_FORMAT") {
            edge_weight_format_ = value;
        } else if (name == "EDGE_DATA_FORMAT") {
            edge_data_format_ = value;
        } else if (name == "NODE_COORD_TYPE") {
            node_coord_type_ = value;
        } else if (name == "DISPLAY_DATA_TYPE") {
            display_data_type_ == value;
        } else {
            is_data_segment = true;
            /* Секция данных */
            FindNextDataSector();
        }

    }

    if (data_sector_type_ == 0) {
        throw std::invalid_argument("Invalid data section");
    }
    return in_file;
}
Matrix *problem_parser::MatrixParseFile(const string& filename) {
    if (!filename.empty()) {
        filename_ = filename;
    } else if (filename_.empty()) {
        throw std::invalid_argument("Filename is not provided");
    }
    if (!started_read) {
        file_point_ = StartParse();
    }
    size_t erase_pos;
    std::string line;

    switch (data_sector_type_) {
        case 1: {
            /* NODE_COORD_SECTION
             * Просто координаты
             * */
            std::string cord;
            size_t i = 0;
            size_t j;
            Matrix cords_matrix = Matrix(dimension_, 2);
            while (std::getline(file_point_, curent_line_)) { // Чтение координат
                line  = curent_line_;
                erase_pos = line.find('\r');
                while (erase_pos != std::string::npos) {
                    line.erase(erase_pos, 1);
                    erase_pos = line.find('\r');
                }
                if (line == "EOF" || line.empty()) {
                    break;
                }
                std::stringstream ss(line);
                SplitLine(ss);
                cords_matrix[i][0] = SplitLine(ss);
                cords_matrix[i][1] = SplitLine(ss);
                i++;
            }
            file_point_.close();

            Matrix *dist_matrix = new Matrix(dimension_);
            for (i = 0; i < dimension_; ++i) {
                (*dist_matrix)[i][i] = 0;
            }

            if (edge_weight_type_ == "EUC_2D") {
                /* XY */
                for (i = 0; i < dimension_ - 1; ++i) {
                    for (j = i + 1; j < dimension_; ++j) {
                        double dist = XYGetDistance(cords_matrix[i][0],
                                                    cords_matrix[i][1],
                                                    cords_matrix[j][0],
                                                    cords_matrix[j][1]);
                        (*dist_matrix)[i][j] = dist;
                        (*dist_matrix)[j][i] = dist;
                    }
                }
            } else if (edge_weight_type_ == "GEO") {
                /* LatLong */
                for (i = 0; i < dimension_ - 1; ++i) {
                    for (j = i + 1; j < dimension_; ++j) {
                        double dist = LatLongGetDistance(cords_matrix[i][0],
                                                         cords_matrix[i][1],
                                                         cords_matrix[j][0],
                                                         cords_matrix[j][1]);
                        (*dist_matrix)[i][j] = dist;
                        (*dist_matrix)[j][i] = dist;
                    }
                }
            } else {
                throw std::invalid_argument("Not supported EDGE_WEIGHT_TYPE");
            }
            return dist_matrix;
        }
        case 8:     {
            /* EDGE_WEIGHT_SECTION
             * Матрица дистанций
             * */
            std::string cord;
            size_t i;
            size_t j;

            Matrix *dist_matrix = new Matrix(dimension_);
            if (edge_weight_format_ == "FULL_MATRIX"){
                /* Полная матрица */
                i = 0;
                j = 0;

                while (std::getline(file_point_, curent_line_)) { // Чтение координат
                    line  = curent_line_;
                    erase_pos = line.find('\r');
                    while (erase_pos != std::string::npos) {
                        line.erase(erase_pos, 1);
                        erase_pos = line.find('\r');
                    }
                    if (line == "EOF" || CheckIfDataSector(curent_line_) > 0) {
                        break;
                    }
                    std::stringstream ss(line);
                    double num = SplitLine(ss);
                    while (!std::isnan(num)){
                        (*dist_matrix)[i][j] = num;
                        num = SplitLine(ss);
                        j++;
                        if (j == dimension_) {
                            i++;
                            j = 0;
                        }
                    }
                }
            } else if (edge_weight_format_ == "UPPER_ROW") {
                /* Верхняя треугольная матрица */
                i = 0;
                j = 1;
                for (int k = 0; k < dimension_; ++k) {
                    (*dist_matrix)[k][k] = 0.0; /* Зануление диагонали */
                }

                while (std::getline(file_point_, curent_line_)) { // Чтение координат
                    line  = curent_line_;
                    erase_pos = line.find('\r');
                    while (erase_pos != std::string::npos) {
                        line.erase(erase_pos, 1);
                        erase_pos = line.find('\r');
                    }
                    if (line == "EOF" || CheckIfDataSector(curent_line_) > 0) {
                        break;
                    }

                    std::stringstream ss(line);
                    double num = SplitLine(ss);
                    while (!std::isnan(num)){
                        (*dist_matrix)[i][j] = num;
                        (*dist_matrix)[j][i] = num;
                        num = SplitLine(ss);
                        j++;
                    }
                    i++;
                    j = i+1;
                }
            } else if (edge_weight_format_ == "LOWER_ROW") {
                /* Нижняя треугольная матрица */
                i = 0;
                j = 1;
                for (int k = 0; k < dimension_; ++k) {
                    (*dist_matrix)[k][k] = 0.0; /* Зануление диагонали */
                }

                while (std::getline(file_point_, curent_line_)) { // Чтение координат
                    line  = curent_line_;
                    erase_pos = line.find('\r');
                    while (erase_pos != std::string::npos) {
                        line.erase(erase_pos, 1);
                        erase_pos = line.find('\r');
                    }
                    if (line == "EOF" || CheckIfDataSector(curent_line_) > 0) {
                        break;
                    }
                    std::stringstream ss(line);
                    double num = SplitLine(ss);;
                    while (!std::isnan(num)){
                        (*dist_matrix)[i][j] = num;
                        (*dist_matrix)[j][i] = num;
                        num = SplitLine(ss);
                        i++;
                    }
                    j++;
                    i = 0;
                }
            } else if (edge_weight_format_ == "UPPER_DIAG_ROW") {
                /* Верхняя треугольная но с диагональю */
                /* TODO: Проверить работоспособность */
                i = -1;
                j = 0;
                while (std::getline(file_point_, curent_line_)) { // Чтение координат
                    line  = curent_line_;
                    erase_pos = line.find('\r');
                    while (erase_pos != std::string::npos) {
                        line.erase(erase_pos, 1);
                        erase_pos = line.find('\r');
                    }
                    if (line == "EOF" || CheckIfDataSector(curent_line_) > 0) {
                        break;
                    }
                    std::stringstream ss(line);
                    double num = SplitLine(ss);;
                    while (!std::isnan(num)){
                        if (num == 0) {
                            i++;
                            j = i;
                            (*dist_matrix)[i][i] = num;
                        } else {
                            j++;
                            (*dist_matrix)[i][j] = num;
                            (*dist_matrix)[j][i] = num;
                        }
                        num = SplitLine(ss);
                    }
                }
            } else if (edge_weight_format_ == "LOWER_DIAG_ROW") {
                /* Нижняя треугольная но с диагональю */
                i = 0;
                j = 0;
                while (std::getline(file_point_, curent_line_)) { // Чтение координат
                    line  = curent_line_;
                    erase_pos = line.find('\r');
                    while (erase_pos != std::string::npos) {
                        line.erase(erase_pos, 1);
                        erase_pos = line.find('\r');
                    }
                    if (line == "EOF" || CheckIfDataSector(curent_line_) > 0) {
                        break;
                    }
                    std::stringstream ss(line);
                    double num = SplitLine(ss);;
                    while (!std::isnan(num)){
                        if (num == 0) {
                            (*dist_matrix)[i][i] = num;
                            i++;
                            j = 0;
                        } else {
                            (*dist_matrix)[i][j] = num;
                            (*dist_matrix)[j][i] = num;
                            j++;
                        }
                        num = SplitLine(ss);
                    }
                }
            }

            file_point_.close();
            return dist_matrix;
        }
        default:
            throw std::invalid_argument("Invalid data sector type");
    }
}

double problem_parser::SplitLine(std::stringstream& ss) {
    std::string cord;
    while (true) {
        // Читаем строку из stringstream
        if (!std::getline(ss, cord, ' ')) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // Проверяем, не пустая ли строка
        if (!cord.empty()) {
            return std::stod(cord);  // Преобразуем строку в число и выводим
        }
    }
}

int problem_parser::CheckIfDataSector(const string& str) {
    if (str == "NODE_COORD_SECTION") {
        return 1;
    } else if (str == "DEPOT_SECTION") {
        return 2;
    } else if (str == "DEMAND_SECTION") {
        return 3;
    } else if (str == "EDGE_DATA_SECTION") {
        return 4;
    } else if (str == "FIXED_EDGES_SECTION") {
        return 5;
    } else if (str == "DISPLAY_DATA_SECTION") {
        return 6;
    } else if (str == "TOUR_SECTION") {
        return 7;
    } else if (str == "EDGE_WEIGHT_SECTION") {
        return 8;
    } else {
        return 0;
    }
}

int problem_parser::FindNextDataSector() {

    string line;
    string delimiter = " : ";
    string delimiter2 = ": ";
    std::string name;
    std::string value;
    size_t pos;
    size_t erase_pos;

    data_sector_type_ == 0;

    while (!curent_line_.empty() && data_sector_type_ == 0) {
        line = curent_line_;
        erase_pos = line.find('\r');
        while (erase_pos != std::string::npos) {
            line.erase(erase_pos, 1);
            erase_pos = line.find('\r');
        }
        erase_pos = line.find(' ');
        while (erase_pos != std::string::npos) {
            line.erase(erase_pos, 1);
            erase_pos = line.find(' ');
        }


        pos = line.find(delimiter);

        if (pos == std::string::npos) {
            pos = line.find(delimiter2);
            name = line.substr(0, pos);
        } else {
            name = line.substr(0, pos);
        }
        data_sector_type_ = CheckIfDataSector(name);
        if (data_sector_type_ == 0){
            std::getline(file_point_, curent_line_);
        }

    }
    return data_sector_type_;
}

