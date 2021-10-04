#ifndef _CRYPTER_H_
#define _CRYPTER_H_
#include <unordered_map>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <codecvt>
#include <sstream>
#include <fstream>
#include <locale>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <regex>
#include <map>

#include "Utils.hpp"

#define _CRYPTER_BEGIN namespace crypter {
#define _CRYPTER_END                     }

#define ZERO 0
#define wcout_t(_Line) std::wcout << _Line << std::endl

#define _STD        ::std::
#define _FILESYSTEM ::std::filesystem::


_CRYPTER_BEGIN
using uint = unsigned int;
using SCS = std::chrono::system_clock;
using pair_wchar = std::pair<wchar_t, uint>;

const uint16_t COUNT_RUS_LETTER = 33;
//? const PATH ------------------------------------------------------ ///
const _FILESYSTEM path KEY = "output/key.txt";

const _FILESYSTEM path ASSAY = "texts/assay.txt";
const _FILESYSTEM path ASSAY_ENCRYPTED = "texts/assay_encrypted.txt";

const _FILESYSTEM path COMPOSITION = "texts/composition.txt";
const _FILESYSTEM path COMPOSITION_ENCRYPTED = "texts/composition_encrypted.txt";

const _FILESYSTEM path RUSSIAN_MONOGRAMS = "sort_rus/russian_monograms.txt";
const _FILESYSTEM path RUSSIAN_BIGRAMS   = "sort_rus/russian_bigrams.txt";
const _FILESYSTEM path RUSSIAN_TRIGRAMS  = "sort_rus/russian_trigrams.txt";

const _FILESYSTEM path RUSSIAN_QUADGRAMS = "rus/russian_quadgrams.txt";
//? const PATH ------------------------------------------------------ ///

inline uint getIndexArray(const uint _Off) { return _Off < 1078 ? _Off - 1072 : _Off == 1105 ? 6 : _Off - 1071; } /// Index: 6 - Russian `E
bool comp(std::pair<std::wstring, double> &_I1, std::pair<std::wstring, double> &_I2) { return _I1.second > _I2.second ; }

class Crypter
{
private:
    uint64_t countAllSymbol;

    std::vector<double> _F1;                                    /// Frequency of monograms
    std::vector<std::vector<double>> _F2;                       /// Frequency of bigrams
    std::vector<std::vector<std::vector<double>>> _F3;          /// Frequency of trigrams

    std::vector<double> _F1_Text;                               /// Frequency of monograms
    std::vector<std::vector<double>> _F2_Text;                  /// Frequency of bigrams
    std::vector<std::vector<std::vector<double>>> _F3_Text;     /// Frequency of trigrams

    void clear();
    void locale() const;
    uint checkLetter(const uint) const;
    uint getLowerID_RU_UTF8(const uint) const;
    _STD wstring getLowerStr(std::wstring) const;

    _STD wstring formatingText(std::wstring &_Text);
    void codeConvertUTF8(std::wfstream &_Stream);

    uint getSeed() const;
    uint getRandom(std::mt19937 &engine, uint, uint) const;
    _STD wstring mutate(std::wstring) const;
    double countDifference();

    void MonoBiTriFileFILL();
    void Normalize(uint64_t _Mono, uint64_t _Bi, uint64_t _Tri);

public:
    Crypter() { locale(); };
    ~Crypter() { };

    _STD wstring EncryptText(std::wstring& key, const std::filesystem::path);
    _STD wstring DecryptText(std::wstring& key, const std::filesystem::path);

    _STD wstring generateRandomKey() const;
    _STD wstring shuffleKey()        const;
    _STD wstring getRussianABC()     const;

    _STD wstring inTextFromFileUTF8 (const std::filesystem::path);
    void outFileUTF8 (_STD wstring &_Text, const std::filesystem::path);

    _STD wstring FrequencyAnalysis(const std::filesystem::path);
};


void Crypter::clear() {
    _F1_Text.clear();
    _F2_Text.clear();
    _F3_Text.clear();

    _F1_Text.resize(COUNT_RUS_LETTER);
    _F2_Text.resize(COUNT_RUS_LETTER, std::vector<double>(COUNT_RUS_LETTER));
    _F3_Text.resize(COUNT_RUS_LETTER, std::vector<std::vector<double>>(COUNT_RUS_LETTER, 
                                                  std::vector<double> (COUNT_RUS_LETTER)));    
}

_STD wstring Crypter::EncryptText(std::wstring& key, const std::filesystem::path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";

    codeConvertUTF8(fin);

    std::wstring text, line;
    while (std::getline(fin, line)) 
        text += line;
    fin.close();

    std::wstring formattedText = formatingText(text);
    std::wstring ABC = getRussianABC();
    std::map<wchar_t, wchar_t> dict;

    for(size_t _Pos = 0; _Pos < ABC.size(); _Pos++)
        dict[ABC[_Pos]] = key[_Pos];

    for(size_t _Pos = 0; _Pos < formattedText.size(); _Pos++)
        formattedText[_Pos] = dict[formattedText[_Pos]];

    return formattedText;
}

_STD wstring Crypter::DecryptText(std::wstring& key, const std::filesystem::path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);
    std::wstring text, line;
    while (std::getline(fin, line)) 
        text += line;
    fin.close();

    std::wstring ABC = getRussianABC();
    std::map<wchar_t, wchar_t> dict;

    for(size_t _Pos = 0; _Pos < ABC.size(); _Pos++)
        dict[key[_Pos]] = ABC[_Pos];

    for(size_t _Pos = 0; _Pos < text.size(); _Pos++)
        text[_Pos] = dict[text[_Pos]];

    return text;
}

_STD wstring Crypter::formatingText(std::wstring &_Text) {
    for (auto &pos : _Text) 
        pos = checkLetter(pos);

    std::wcmatch charRg;    
    std::wregex rg (L"[а-яё]");
    std::wostringstream wostream;
    while (std::regex_search(_Text.c_str(), charRg, rg)) {
        wostream << charRg.str();
        _Text = charRg.suffix().str();
    }
    return wostream.str();
}

_STD wstring Crypter::generateRandomKey() const {
    std::mt19937 engine(getSeed());

    std::wstring ABC = getRussianABC();
    for(size_t _Pos = ABC.size() - 1; _Pos > 0; _Pos--) {
        uint rdPos = getRandom(engine, ZERO, _Pos + 1);
        std::swap(ABC[_Pos], ABC[rdPos]);
    }
    return ABC;
}

_STD wstring Crypter::shuffleKey() const {
    std::wstring ABC = getRussianABC();
    std::random_device rd;
    std::mt19937 engine(rd());
 
    std::shuffle(ABC.begin(), ABC.end(), engine);
    return ABC;
}

_STD wstring Crypter::getRussianABC() const {
    const size_t LET_A = 1072;
    const size_t POS_E = 6, ID = 1105;

    std::wstring ABC;
    ABC.resize(COUNT_RUS_LETTER);

    for (size_t _Pos = 0; _Pos < POS_E; _Pos++) 
        ABC[_Pos] = LET_A + _Pos;
    ABC[POS_E] = ID;
    for (size_t _Pos = POS_E; _Pos < COUNT_RUS_LETTER - 1; _Pos++) 
        ABC[_Pos + 1] = LET_A + _Pos;

    return ABC;
}

uint Crypter::getRandom(std::mt19937 &engine, uint _minValue, uint _maxValue) const {
    return (engine() % (_maxValue - _minValue) + _minValue); 
}
uint Crypter::checkLetter(const uint _Id) const {
    return ((_Id >= 1072 && _Id <= 1103) || _Id == 1105 ? _Id : getLowerID_RU_UTF8(_Id));
} 
uint Crypter::getLowerID_RU_UTF8(const uint _Id) const { return _Id == 1025 ? 1105 : _Id + 32; }
uint Crypter::getSeed() const { return SCS::now().time_since_epoch().count(); }
void Crypter::locale() const { setlocale(LC_ALL, "Russian"); }  

_STD wstring Crypter::getLowerStr(_STD wstring _Str) const {
    for(auto &_Symbol : _Str)
        _Symbol = checkLetter(_Symbol);
    return _Str;
}


_STD wstring Crypter::inTextFromFileUTF8(const std::filesystem::path _Path) {
    std::wfstream _Fin(_Path, std::wfstream::in); 
    if(!_Fin) 
        std::cerr << "FILE WAS NOT OPEN!";

    codeConvertUTF8(_Fin);
    std::wstring _Text, _Line;
    std::wstringstream _Stream;
    while (std::getline(_Fin, _Line)) 
        _Stream << _Line;
    _Text = _Stream.str();
    _Fin.close();

    return _Text;
}

void Crypter::outFileUTF8(std::wstring& _Text, const std::filesystem::path _Path) {
    std::wfstream fout(_Path, std::wfstream::out);
    codeConvertUTF8(fout);
    fout << _Text;
    fout.close();
}

void Crypter::codeConvertUTF8(std::wfstream &_Stream) { 
    _Stream.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>)); 
}

_STD wstring Crypter::mutate(std::wstring _Key) const {
    std::wstring _K = _Key;
    std::random_device rd;
    std::mt19937 engine(rd());

    std::swap(_K[getRandom(engine, 0, _K.size())], 
              _K[getRandom(engine, 0, _K.size())]);
    return _K;
}

void Crypter::MonoBiTriFileFILL() {
    _F1.resize(COUNT_RUS_LETTER);
    _F2.resize(COUNT_RUS_LETTER, std::vector<double>(COUNT_RUS_LETTER));
    _F3.resize(COUNT_RUS_LETTER, std::vector<std::vector<double>>(COUNT_RUS_LETTER, 
                                             std::vector<double> (COUNT_RUS_LETTER)));

    std::wstring _Str, _Line;
    uint64_t _Count;

    uint64_t countMonograms = 0, countBigrams = 0, countTrigrams = 0; 

    std::wfstream _Fin(RUSSIAN_MONOGRAMS, std::wfstream::in); 
    codeConvertUTF8(_Fin);

    while (std::getline(_Fin, _Line)) {
        std::wstringstream _Stream(_Line);
        _Stream >> _Str >> _Count;

        countMonograms += _Count;
        std::wstring lower = getLowerStr(_Str);

        _F1 [getIndexArray(lower[0])] = _Count;
    }
    _Fin.close();
    

    _Fin.open(RUSSIAN_BIGRAMS, std::wfstream::in); 
    codeConvertUTF8(_Fin);

    while (std::getline(_Fin, _Line)) {
        std::wstringstream _Stream(_Line);
        _Stream >> _Str >> _Count;

        countBigrams += _Count;
        std::wstring lower = getLowerStr(_Str);

        _F2 [getIndexArray(lower[0])]
            [getIndexArray(lower[1])] = _Count;
    }
    _Fin.close();


    _Fin.open(RUSSIAN_TRIGRAMS, std::wfstream::in); 
    codeConvertUTF8(_Fin);

    while (std::getline(_Fin, _Line)) {
        std::wstringstream _Stream(_Line);
        _Stream >> _Str >> _Count;

        countTrigrams += _Count;
        std::wstring lower = getLowerStr(_Str);

        _F3 [getIndexArray(lower[0])]
            [getIndexArray(lower[1])]
            [getIndexArray(lower[2])] = _Count;
    }
    _Fin.close();

    Normalize(countMonograms, countBigrams, countTrigrams);
}

void Crypter::Normalize(uint64_t _Mono, uint64_t _Bi, uint64_t _Tri) {
    const size_t _Count = COUNT_RUS_LETTER;

    for (size_t i = 0; i < _Count; i++)
        _F1[i] /= _Mono;

    for (size_t i = 0; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++)
            _F2[i][j] /= _Bi;

    for (size_t i = 0; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++)
            for (size_t k = 0; k < _Count; k++)
                _F3[i][j][k] /= _Tri;
}

_STD wstring Crypter::FrequencyAnalysis(const std::filesystem::path _Path) {
    const size_t _Count = COUNT_RUS_LETTER;

    clear();                                                                                    /// Очищаем вектора
    std::wstring text = inTextFromFileUTF8(_Path);                                              /// Получаем текст, который нужно расшифровать
    MonoBiTriFileFILL();                                                                        /// Заполняем Моно-Би-Триграмы из файлов

/// Monograms
    std::vector<std::pair<wchar_t, uint>> countMonograms(_Count);
    std::wstring ABC = getRussianABC();
    for (size_t _Pos = 0; _Pos < ABC.size(); _Pos++) 
        countMonograms[_Pos] = std::pair(ABC[_Pos], ZERO);
    
    for (const auto &_El : text) 
        countMonograms[getIndexArray(_El)].second++;
    
    std::sort(countMonograms.begin(), countMonograms.end(), 
        [] (const std::pair<wchar_t, uint> &_Left, 
            const std::pair<wchar_t, uint> &_Right) { return _Left.second > _Right.second; });  /// Сортировка в порядке убывания

    countAllSymbol = std::accumulate(countMonograms.begin(), countMonograms.end(), 0, 
        [](const size_t _S, const auto &_Elem) { return _S + _Elem.second; });                  /// Количество всех символов текста

    std::wstring monogramsKey;                                                                  /// Ключ отсортированный в порядке 
    std::wstringstream w_stream;                                                                /// количества наибольших вхождений
    for (auto &_El : countMonograms)                                                            /// к наименьшим
        w_stream << _El.first;
    monogramsKey = w_stream.str();

    std::wfstream fin(RUSSIAN_MONOGRAMS, std::wfstream::in);                                    /// Открываем файл с частотой вхождения русских букв
    codeConvertUTF8(fin);                                                                       /// Устанавливаем кодировку UTF-8 
    
    std::vector<std::pair<wchar_t, uint>> countMonogramsFile(_Count);                           /// Количество частоты вхождений букв из файла
    std::wstring _Str, line;
    for (size_t _Pos = 0, _Size; std::getline(fin, line); _Pos++) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;

        countMonogramsFile[_Pos] = std::pair(getLowerStr(_Str)[0], _Size);
    }
    fin.close();

    std::sort(countMonogramsFile.begin(), countMonogramsFile.end(),                             /// Сортируем вектор, для дальнейшего
        [] (const std::pair<wchar_t, uint> &_Left,                                              /// сопоставления пары wchar_t : wchar_t
            const std::pair<wchar_t, uint> &_Right) { return _Left.second > _Right.second; }); 


    std::vector<std::pair<wchar_t, wchar_t>> monogramsPair(_Count);                             /// Вектор : Ключ <-> Значение
    for (size_t _Pos = 0; _Pos < monogramsKey.size(); _Pos++) 
        monogramsPair[getIndexArray(monogramsKey[_Pos])] = std::pair(monogramsKey[_Pos], countMonogramsFile[_Pos].first);
    

    w_stream = std::wstringstream();
    for (const auto &x : text) 
        w_stream << monogramsPair[getIndexArray(x)].second;                                     /// Заменяем каждую букву текста на  
    std::wstring textDecrypted = w_stream.str();                                                /// соответсвующую ей частотности в векторе


/// Bigrams
    for (int _Pos = 0; _Pos < textDecrypted.size() - 1; _Pos++) {
        std::wstring currentBigram = textDecrypted.substr(_Pos, 2);

        _F2_Text[getIndexArray(currentBigram[0])]
                [getIndexArray(currentBigram[1])]++;
    }

/// Trigrams
    for (int _Pos = 0; _Pos < textDecrypted.size() - 2; _Pos++) {
        std::wstring currentTrigram = textDecrypted.substr(_Pos, 3);

        _F3_Text[getIndexArray(currentTrigram[0])]
                [getIndexArray(currentTrigram[1])]
                [getIndexArray(currentTrigram[2])]++;
    }

    for (const auto &_El : countMonograms)                                                      /// Помещаем частоту вхождения каждой 
        _F1_Text[getIndexArray(_El.first)] = static_cast<double>(_El.second) / countAllSymbol;  /// буквы в отдельный вектор

    for (size_t i = 0, countBigrams = countAllSymbol - 1; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++) 
            _F2_Text[i][j] /= countBigrams;

    for (size_t i = 0, countTrigrams = countAllSymbol - 2; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++) 
                for (size_t k = 0; k < _Count; k++) 
                    _F3_Text[i][j][k] /= countTrigrams;


/// ----------------------------------------------------------------------- ///
/// CONST GENETIC ALGORITHM
    const uint populationSize = 50;             /// Количество популяций
    const uint bestPopulation = 10;             /// Количество отбираемых особей
    const uint generationCount = 100;           /// Количество генераций поколений
/// CONST GENETIC ALGORITHM

    double _Efficiency = countDifference();

    std::vector<std::pair<std::wstring, double>> populationEfficiency(populationSize);          /// Популяция (КЛЮЧ - ЕФФЕКТИВНОСТЬ)
    for (size_t i = 0; i < populationSize; i++) 
        populationEfficiency[i] = std::pair(monogramsKey, _Efficiency);        

    std::vector<std::pair<std::wstring, double>> tempPopulationEfficiency;  


    for (int _Generation = 0; _Generation < generationCount; _Generation++) {
        std::wcout << "\r"<< _Generation + 1 << L"-ая итерация" << "\t" << "Efficiency: " << populationEfficiency[0].second;

        std::sort(populationEfficiency.begin(), populationEfficiency.end(), 
            [] (const std::pair<std::wstring, double> &_Left, 
                const std::pair<std::wstring, double> &_Right) { return _Left.second < _Right.second; });

        tempPopulationEfficiency.clear();
        for (int i = 0; i < bestPopulation; i++) 
            tempPopulationEfficiency.push_back(populationEfficiency[i]);

        int descendantCount = (populationSize - bestPopulation) / bestPopulation;
        for (int i = 0; i < bestPopulation; i++) 
            for (int j = 0; j < descendantCount; j++) {
                std::wstring _Key = mutate(populationEfficiency[i].first);
                tempPopulationEfficiency.push_back(std::pair(_Key, ZERO));
            }
        

        for (int i = bestPopulation; i < tempPopulationEfficiency.size(); i++) {
            monogramsKey = tempPopulationEfficiency[i].first;

            monogramsPair.clear();
            for (size_t _Pos = 0; _Pos < monogramsKey.size(); _Pos++) 
                monogramsPair[getIndexArray(monogramsKey[_Pos])] = std::pair(monogramsKey[_Pos], countMonogramsFile[_Pos].first);


            w_stream = std::wstringstream();
            for (const auto &x : text) 
                w_stream << monogramsPair[getIndexArray(x)].second; 
            textDecrypted = w_stream.str(); 

            clear();                                                                                        /// Очистка векторов
            std::wstring _Сurrent;
            uint _Index1, _Index2, _Index3;
            for (int _Pos = 0; _Pos < textDecrypted.size() - 2; _Pos++) {
                _Сurrent = textDecrypted.substr(_Pos, 3);

                _Index1 = getIndexArray(_Сurrent[0]);
                _Index2 = getIndexArray(_Сurrent[1]);
                _Index3 = getIndexArray(_Сurrent[2]);

                _F1_Text[_Index1]++;                                                                        /// MONOGRAMS
                _F2_Text[_Index1][_Index2]++;                                                               /// BIGRAMS
                _F3_Text[_Index1][_Index2][_Index3]++;                                                      /// TRIGRAMS
            }
            _F1_Text[_Index2]++;    
            _F1_Text[_Index3]++;
            _F2_Text[_Index2][_Index3]++;


            for (size_t i = 0; i < _Count; i++)
                    _F1_Text[i] /= countAllSymbol;

            for (size_t i = 0, countBigrams = countAllSymbol - 1; i < _Count; i++)
                for (size_t j = 0; j < _Count; j++) 
                    _F2_Text[i][j] /= countBigrams;

            for (size_t i = 0, countTrigrams = countAllSymbol - 2; i < _Count; i++)
                for (size_t j = 0; j < _Count; j++) 
                        for (size_t k = 0; k < _Count; k++) 
                            _F3_Text[i][j][k] /= countTrigrams;

            _Efficiency = countDifference();

            tempPopulationEfficiency[i].second = _Efficiency;
        }
        populationEfficiency = tempPopulationEfficiency;

    }

    std::sort(populationEfficiency.begin(), populationEfficiency.end(),                         /// Сортируем вектор по эффективности
        [] (const std::pair<std::wstring, double> &_Left, 
            const std::pair<std::wstring, double> &_Right) { return _Left.second < _Right.second; });
    monogramsKey = populationEfficiency[0].first;

    std::wstring bestKey;
    bestKey.resize(_Count);

    size_t _Pos = 0;
    for (size_t _Pos = 0; _Pos < _Count; _Pos++) 
        bestKey[getIndexArray(countMonogramsFile[_Pos].first)] = monogramsKey[_Pos];

    return bestKey;
}


double metric(double _F, double _T, int64_t& _Count) {
    if (_F == 0 && _T != 0) {
        _Count++; 
        return 1;
    } else if (_T != 0) {
        _Count++;
        return fabs(_F - _T) / std::max(_F, _T);
    } else return 0;
}

double Crypter::countDifference() {
    const size_t _Count = COUNT_RUS_LETTER;
    int64_t count1 = 0, count2 = 0, count3 = 0;
    double result1 = 0, result2 = 0, result3 = 0;

    for (size_t i = 0; i < _Count; i++)
        result1 += metric(_F1[i], _F1_Text[i], count1);

    for (size_t i = 0; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++)
            result2 += metric(_F2[i][j], _F2_Text[i][j], count2);

    for (size_t i = 0; i < _Count; i++)
        for (size_t j = 0; j < _Count; j++)
            for (size_t k = 0; k < _Count; k++)
                result3 += metric(_F3[i][j][k], _F3_Text[i][j][k], count3);

    result1 /= count1;
    result2 /= count2;
    result3 /= count3;

    double weight1 = 1, weight2 = 2, weight3 = 3;           /// Веса 
    double weightSum = weight1 + weight2 + weight3;         /// Самыми приоритетными являются Триграммы 

    return ( 
        weight1 * result1 + 
        weight2 * result2 + 
        weight3 * result3 
    ) / weightSum;
}
_CRYPTER_END
#undef _STD
#endif // _CRYPTER_H_