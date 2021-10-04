#ifndef _UTILS_H_
#define _UTILS_H_
#include <iostream>
#include <locale>
#include <fstream>
#include <codecvt>
#include <vector>
#include <sstream>
#include <string>


#define _UTILS_BEGIN namespace utils {
#define _UTILS_END                   }
const int COUNT_RUS_LET = 33;

_UTILS_BEGIN

const std::string RUSSIAN_MONOGRAMS = "rus/russian_monograms.txt";
const std::string RUSSIAN_BIGRAMS   = "rus/russian_bigrams.txt";
const std::string RUSSIAN_TRIGRAMS  = "rus/russian_trigrams.txt";

int getIndex(const size_t _Off) { return _Off == 1025 ? 1046 : _Off < 1046 ? _Off : _Off + 1 ; } /// Index: 1025 - Russian `E
int getLetter(const size_t _Off) { return _Off == 6 ? 1025 : _Off < 6 ? _Off + 1040 : _Off + 1039 ; } /// Index: 1025 - Russian `E



int getLetterLower(const size_t _Off) { return _Off == 6 ? 1105 : _Off < 6 ? _Off + 1072 : _Off + 1071 ; } /// Index: 1025 - Russian `E
wchar_t getLowerID_RU_UTF8(const int _Id) { return _Id == 1025 ? 1105 : _Id + 32; }


std::wstring getStrBI(int i, int j) {

    wchar_t letter1 = getLetter(i);
    wchar_t letter2 = getLetter(j);

    std::wstringstream w_stream;
    w_stream << letter1 << letter2;
    std::wstring _Str(w_stream.str());

    return _Str;
}

std::wstring getStrBILower(int i, int j) {

    wchar_t letter1 = getLetterLower(i);
    wchar_t letter2 = getLetterLower(j);

    std::wstringstream w_stream;
    w_stream << letter1 << letter2;
    std::wstring _Str(w_stream.str());

    return _Str;
}

std::wstring getStrTRI(int i, int j, int k) {

    wchar_t letter1 = getLetter(i);
    wchar_t letter2 = getLetter(j);
    wchar_t letter3 = getLetter(k);

    std::wstringstream w_stream;
    w_stream << letter1 << letter2 << letter3;
    std::wstring _Str(w_stream.str());

    return _Str;
}

std::wstring getStrTRILower(int i, int j, int k) {

    wchar_t letter1 = getLetterLower(i);
    wchar_t letter2 = getLetterLower(j);
    wchar_t letter3 = getLetterLower(k);

    std::wstringstream w_stream;
    w_stream << letter1 << letter2 << letter3;
    std::wstring _Str(w_stream.str());

    return _Str;
}

class Utils
{
private:
    void locale() const;
    void codeConvertUTF8(std::wfstream &_Stream);


public:
     Utils() { locale(); };
    ~Utils() { };

    void MONO();
    void BI();
    void TRI();

};

void Utils::codeConvertUTF8(std::wfstream &_Stream) { _Stream.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>)); }
void Utils::locale() const { setlocale(LC_ALL, "Russian"); }  

void Utils::MONO() {
    std::wfstream fin(RUSSIAN_MONOGRAMS, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);

    std::vector<std::pair<std::wstring, uint64_t>> _monograms;
    _monograms.resize(COUNT_RUS_LET);

    std::wstring _Str, line; 
    uint64_t _Size;
    size_t _Pos = 0;

    for (size_t _Pos = 0; std::getline(fin, line); _Pos++) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;
        _monograms[_Pos] = std::pair(_Str, _Size);
    }

    std::sort(_monograms.begin(), _monograms.end(), 
        [](const std::pair<std::wstring, uint64_t> &_Left, const std::pair<std::wstring, uint64_t> &_Right) 
            { return getIndex(_Left.first.at(0)) < getIndex(_Right.first.at(0)); });

    std::wfstream fout("sort_lower_rus/russian_monograms.txt", std::wfstream::out);
    codeConvertUTF8(fout);

    std::wstringstream w_stream;
    for (const auto _Elem : _monograms) 
        w_stream << getLowerID_RU_UTF8(_Elem.first.at(0)) << " " << _Elem.second << std::endl;      /// Тут убрать getLowerID_RU_UTF8 если нужны большие

    fout << w_stream.str();
    fout.close();
}

void Utils::BI() {
    std::wfstream fin(RUSSIAN_BIGRAMS, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);

    std::vector<std::pair<std::wstring, uint64_t>> _bigrams;
    std::wstring _Str, line; 
    uint64_t _Size;
    size_t _Pos = 0;

    for (size_t _Pos = 0; std::getline(fin, line); _Pos++) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;

        _bigrams.push_back(std::pair(_Str, _Size));
    }



    std::vector<std::pair<std::wstring, uint64_t>> _sortedBigrams;
    _sortedBigrams.resize(COUNT_RUS_LET * COUNT_RUS_LET);

    for (int i = 0; i < COUNT_RUS_LET; i++) {
        for (int j = 0; j < COUNT_RUS_LET; j++) {
            _Pos = i * COUNT_RUS_LET + j;

            std::wstring _Bi = getStrBI(i, j);
            auto got = find_if(_bigrams.cbegin(), _bigrams.cend(), 
                [&_Bi](const std::pair<std::wstring, uint64_t> &_I){ return _I.first == _Bi; });

            if(got == _bigrams.end()) 
                _sortedBigrams[_Pos] = std::pair(getStrBILower(i, j), 0);               /// Тут менять если на большие нужно
            else 
                _sortedBigrams[_Pos] = std::pair(getStrBILower(i, j), got->second);     /// Тут менять если на большие нужно
        }
    }

    std::wfstream fout("sort_lower_rus/russian_bigrams.txt", std::wfstream::out);
    codeConvertUTF8(fout);

    std::wstringstream w_stream;
    for (const auto _Elem : _sortedBigrams) 
        w_stream << _Elem.first << " " << _Elem.second << std::endl;

    fout << w_stream.str();
    fout.close();
}


void Utils::TRI() {
    std::wfstream fin(RUSSIAN_TRIGRAMS, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);

    std::vector<std::pair<std::wstring, uint64_t>> _trigrams;
    std::wstring _Str, line; 
    uint64_t _Size;
    size_t _Pos = 0;

    for (size_t _Pos = 0; std::getline(fin, line); _Pos++) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;
        _trigrams.push_back(std::pair(_Str, _Size));
    }

    std::vector<std::pair<std::wstring, uint64_t>> _sortedTrigrams;
    _sortedTrigrams.resize(COUNT_RUS_LET * COUNT_RUS_LET * COUNT_RUS_LET);

    for (int i = 0; i < COUNT_RUS_LET; i++) {
        for (int j = 0; j < COUNT_RUS_LET; j++) {
            for (int k = 0; k < COUNT_RUS_LET; k++) {

                _Pos = ( i * COUNT_RUS_LET + j ) * COUNT_RUS_LET + k;

                std::wstring _Tri = getStrTRI(i, j, k);
                auto got = find_if(_trigrams.cbegin(), _trigrams.cend(), 
                    [&_Tri](const std::pair<std::wstring, uint64_t> &_I){ return _I.first == _Tri; });

                if(got == _trigrams.end()) 
                    _sortedTrigrams[_Pos] = std::pair(getStrTRILower(i, j, k), 0);              /// Тут менять если на большие нужно
                else 
                    _sortedTrigrams[_Pos] = std::pair(getStrTRILower(i, j, k), got->second);    /// Тут менять если на большие нужно
            }
        }
    }

    std::wfstream fout("sort_lower_rus/russian_trigrams.txt", std::wfstream::out);
    codeConvertUTF8(fout);

    std::wstringstream w_stream;
    for (const auto _Elem : _sortedTrigrams) 
        w_stream << _Elem.first << " " << _Elem.second << std::endl;

    fout << w_stream.str();
    fout.close();
}
_UTILS_END

namespace Timer {

    using std::chrono::high_resolution_clock;
    using std::chrono::time_point;
    using std::chrono::duration;

    class Timer
    {
    private:
        time_point<high_resolution_clock> timeStart;
        time_point<high_resolution_clock> timeEnd;
    public:
        Timer() : timeStart(high_resolution_clock::now()){ }
        ~Timer() { };

        void setTimeStart() { timeStart = high_resolution_clock::now(); };

        void setTimeEnd();
        float getElapsed() const;

        friend std::ostream& operator<< (std::ostream& out, Timer& point) {
            return out << "SECONDS: " << point.getElapsed() << std::endl;
        }
    };
    void  Timer::setTimeEnd()       { timeEnd = high_resolution_clock::now(); }
    float Timer::getElapsed() const { return duration<float>(timeEnd - timeStart).count(); }

}

namespace ProgressBar {
    class ProgressBarConsole
    {
    private:
        std::string _progressBar;
        uint8_t     _minRange;
        uint8_t     _Position;
        uint8_t     _maxRange;
        const char  symbolRun = '*';    


    public:
        ProgressBarConsole(const uint8_t maxRange = 10);
        void setRange(const uint8_t maxRange);
        void stepIt(uint8_t _Cur);

    private:
        void init();
        void UpdateProgress();
    };

    ProgressBarConsole::ProgressBarConsole(const uint8_t maxRange) 
        : _minRange(0), _maxRange(maxRange) { init(); }

    void ProgressBarConsole::init() {
        _Position = 1;

        _progressBar.clear();
        _progressBar.resize(20 + 2, static_cast<char>(32));
        _progressBar[_minRange] = '[';
        _progressBar[20 + 1] = ']';

        UpdateProgress();
    }

    void ProgressBarConsole::stepIt(uint8_t _Cur) {
        float _Fr = static_cast<float>(_Cur + 1) / _maxRange;
        std::string _Str = std::to_string(_Fr);

        const size_t _Count = 5;
        for (size_t _Pos = 0; _Pos < _Count; _Pos++)
            _progressBar[_Pos + 1] = _Str[_Pos];

        UpdateProgress();
    }

    void ProgressBarConsole::UpdateProgress() { std::cout << _progressBar << "\r"; }

    void ProgressBarConsole::setRange(const uint8_t maxRange) { 
        _maxRange = maxRange; 
        init();
    }
}


#endif // _UTILS_H_