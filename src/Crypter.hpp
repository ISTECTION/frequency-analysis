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

const _FILESYSTEM path RUSSIAN_MONOGRAMS = "rus/russian_monograms.txt";
const _FILESYSTEM path RUSSIAN_BIGRAMS = "rus/russian_bigrams.txt";
const _FILESYSTEM path RUSSIAN_TRIGRAMS = "rus/russian_trigrams.txt";
const _FILESYSTEM path RUSSIAN_QUADGRAMS = "rus/russian_quadgrams.txt";
//? const PATH ------------------------------------------------------ ///

uint getIndexArray(const uint _Off) { return _Off < 1078 ? _Off - 1072 : _Off == 1105 ? 6 : _Off - 1071; } /// Index: 6 - Russian `E
bool comp(std::pair<std::wstring, double> &_I1, std::pair<std::wstring, double> &_I2) { return _I1.second > _I2.second ; }

class Crypter
{
private:
    uint countAllSymbol;

    void locale() const;
    uint checkLetter(const uint) const;
    uint getLowerID_RU_UTF8(const uint) const;
    _STD wstring getLowerStr(std::wstring) const;

    _STD wstring formatingText(std::wstring &_Text);
    void codeConvertUTF8(std::wfstream &_Stream);

    uint getSeed() const;
    uint getRandom(std::mt19937 &engine, uint, uint) const;
    _STD wstring mutate(std::wstring) const;
    double countDifference(std::unordered_map<std::wstring, double> &monogramsText,
                           std::unordered_map<std::wstring, double> &monogramsFile,
                           std::unordered_map<std::wstring, double> &bigramsText,
                           std::unordered_map<std::wstring, double> &bigramsFile,
                           std::unordered_map<std::wstring, double> &trigramsText,
                           std::unordered_map<std::wstring, double> &trigramsFile);

    std::unordered_map<std::wstring, double>& inFileUTF8_FORMAT_WCHAR_COUNT(std::unordered_map<std::wstring, double> &, const std::filesystem::path);
public:
    Crypter() { locale(); };
    ~Crypter() { };

    _STD wstring EncryptText(std::wstring& key, const std::filesystem::path);
    _STD wstring DecryptText(std::wstring& key, const std::filesystem::path);

    _STD wstring generateRandomKey() const;
    _STD wstring shuffleKey()        const;
    _STD wstring getRussianABC()     const;

    _STD wstring inKeyFromFileUTF8 (const std::filesystem::path);
    void outFileUTF8 (_STD wstring &_Text, const std::filesystem::path);

    void FrequencyAnalysis(const std::filesystem::path);
};

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


std::unordered_map<std::wstring, double>& Crypter::inFileUTF8_FORMAT_WCHAR_COUNT(std::unordered_map<std::wstring, double> &_Map,
    const std::filesystem::path _Path) {
    std::wfstream fin(_Path, std::wfstream::in);
    if (!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);

    std::wstring _Str, line; 
    uint _Size;
    while (std::getline(fin, line)) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;
        _Map[getLowerStr(_Str)] = _Size;
    }
    fin.close();

    unsigned long long countAll__Grams = std::accumulate(_Map.begin(), _Map.end(), 0, 
        [](const uint _S, const auto &_Elem) { return _S + _Elem.second; });
    
    for (auto &_Elem : _Map) 
        _Elem.second /= countAll__Grams;
    
    return _Map;
}

_STD wstring Crypter::inKeyFromFileUTF8(const std::filesystem::path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if (!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);
    std::wstring key;
    std::getline(fin, key);
    fin.close();
    return key;
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
    std::random_device rd;
    std::mt19937 engine(rd());
    std::wstring _K = _Key;

    std::swap(_K[getRandom(engine, 0, _K.size())], 
              _K[getRandom(engine, 0, _K.size())]);
    return _K;
}


void Crypter::FrequencyAnalysis(const std::filesystem::path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";

    codeConvertUTF8(fin);
    std::wstring text, line;
    std::wstringstream w_stream;
    while (std::getline(fin, line)) 
        w_stream << line;
    text = w_stream.str();
    fin.close();


/// Monograms
    std::vector<std::pair<wchar_t, uint>> countMonograms;
    countMonograms.resize(COUNT_RUS_LETTER);

    std::wstring ABC = getRussianABC();
    for (size_t _Pos = 0; _Pos < ABC.size(); _Pos++) 
        countMonograms[_Pos] = std::pair(ABC[_Pos], ZERO);
    
    for (const auto &_El : text) 
        countMonograms[getIndexArray(_El)].second++;
    
    std::sort(countMonograms.begin(), countMonograms.end(), 
        [] (const pair_wchar &_Left, const pair_wchar &_Right) { return _Left.second > _Right.second; });

    countAllSymbol = std::accumulate(countMonograms.begin(), countMonograms.end(), 0, 
        [](const size_t _S, const auto &_Elem) { return _S + _Elem.second; });


    /// TODO:
        std::unordered_map<std::wstring, double> monogramsText;
        for (auto &x : countMonograms) {
            std::wstring temp(1, x.first);
            monogramsText[temp] = static_cast<double>(x.second) / countAllSymbol;
        }
    /// TODO:


    w_stream = std::wstringstream();
    std::wstring monogramsKey;
    for(auto &_El : countMonograms)
        w_stream << _El.first;
    monogramsKey = w_stream.str();


    std::unordered_map<wchar_t, wchar_t> monogramsMap;
    fin.open(RUSSIAN_MONOGRAMS, std::wfstream::in);
    codeConvertUTF8(fin);

    /// TODO:
        std::unordered_map<std::wstring, double> monogramsFileMap;
    /// TODO:

    std::wstring fileKey;
    std::wstring _Str;
    uint _Size;
    for (size_t _Pos = 0; std::getline(fin, line); _Pos++) {
        std::wstringstream w_stream(line);
        w_stream >> _Str >> _Size;

        /// TODO:
            monogramsFileMap[getLowerStr(_Str)] = _Size;
        /// TODO:

        monogramsMap[monogramsKey[_Pos]] = getLowerStr(_Str).at(0);
        fileKey += getLowerStr(_Str).at(0);
    }
    fin.close();


    /// TODO:
        unsigned long long countMonogramsFileMap = std::accumulate(monogramsFileMap.begin(), monogramsFileMap.end(), 0, 
            [](const uint _S, const auto &_Elem) { return _S + _Elem.second; });
        
        for (auto &_Elem : monogramsFileMap) 
            _Elem.second /= countMonogramsFileMap;
    /// TODO:

    w_stream = std::wstringstream();
    for (auto &x : text) 
        w_stream << monogramsMap[x];
    std::wstring textDecrypted = w_stream.str();

/// Reserve key
    std::wstring reservMonogramsKey = monogramsKey;

/// Bigrams
    std::unordered_map<std::wstring, double> bigramsText;
    for (int _Pos = 0; _Pos < textDecrypted.size() - 1; _Pos++) {
        std::wstring currentBigram = textDecrypted.substr(_Pos, 2);

        if(bigramsText.find(currentBigram) == bigramsText.end()) 
            bigramsText[currentBigram] = 1;
        else
            bigramsText[currentBigram]++;
    }
    int countBigrams = countAllSymbol - 1;
    for (auto &_Elem : bigramsText)
        _Elem.second /= countBigrams;

    std::unordered_map<std::wstring, double> bigramsFileMap;
    bigramsFileMap = inFileUTF8_FORMAT_WCHAR_COUNT(bigramsFileMap, RUSSIAN_BIGRAMS);

/// Trigrams
    std::unordered_map<std::wstring, double> trigramsText;
    for (int _Pos = 0; _Pos < textDecrypted.size() - 2; _Pos++) {
        std::wstring currentTrigram = textDecrypted.substr(_Pos, 3);

        if(trigramsText.find(currentTrigram) == trigramsText.end()) 
            trigramsText[currentTrigram] = 1;
        else
            trigramsText[currentTrigram]++;
    }
    int countTrigrams = countAllSymbol - 2;
    for (auto &_Elem : trigramsText)
        _Elem.second /= countTrigrams;

    std::unordered_map<std::wstring, double> trigramsFileMap;
    trigramsFileMap = inFileUTF8_FORMAT_WCHAR_COUNT(trigramsFileMap, RUSSIAN_TRIGRAMS);


/// ----------------------------------------------------------------------- ///
/// CONST GENETIC ALGORITHM
    const uint populationSize = 50;      /// Количество популяций
    const uint bestPopulation = 10;      /// Количество отбираемых особей
    const uint generationCount = 300;    /// Количество генераций поколений
/// CONST GENETIC ALGORITHM

    double _Efficiency = countDifference(monogramsText, monogramsFileMap, bigramsText, bigramsFileMap, trigramsText, trigramsFileMap);


    std::vector<std::pair<std::wstring, double>> populationEfficiency;
    populationEfficiency.resize(populationSize);
    for (size_t i = 0; i < populationSize; i++) 
        populationEfficiency[i] = std::pair(monogramsKey, _Efficiency);
        
    std::vector<std::pair<std::wstring, double>> tempPopulationEfficiency;

    for (int _Generation = 0; _Generation < 100; _Generation++) {
        std::wcout << L"Я на " << _Generation << L"-ой итерации!" << std::endl;

        std::sort(populationEfficiency.begin(), populationEfficiency.end(), 
            [] (const std::pair<std::wstring, double> &_Left, 
                const std::pair<std::wstring, double> &_Right) { return _Left.second < _Right.second; });


        tempPopulationEfficiency.clear();
        for (int i = 0; i < bestPopulation; i++) 
            tempPopulationEfficiency.push_back(populationEfficiency[i]);

        int descendantCount = (populationSize - bestPopulation) / bestPopulation;
        for (int i = 0; i < bestPopulation; i++) 
            for (int j = 0; j < descendantCount; j++) 
            {
                std::wstring _Key = mutate(populationEfficiency[i].first);
                tempPopulationEfficiency.push_back(std::pair(_Key, ZERO));
            }
        
        int descendantResidue = populationSize - bestPopulation - descendantCount * bestPopulation;
        for (int i = 0; i < descendantResidue; i++) {
            std::wstring _Key = mutate(populationEfficiency[i].first);
            tempPopulationEfficiency.push_back(std::pair(_Key, ZERO));
        }

        for (int i = 0; i < tempPopulationEfficiency.size(); i++) {
            std::wcout << populationEfficiency[i].first << " " << populationEfficiency[i].second << std::endl;

            monogramsKey = tempPopulationEfficiency[i].first;

            monogramsMap.clear();
            int _Pos = 0;
            for (const auto &x: monogramsKey) {
                monogramsMap[x] = fileKey[_Pos];
                _Pos++;
            }
                

            w_stream = std::wstringstream();
            for (auto &x : text) 
                w_stream << monogramsMap[x];
            textDecrypted = w_stream.str();


            monogramsText.clear();
            for (int _Pos = 0; _Pos < textDecrypted.size(); _Pos++) {
                std::wstring currentBigram = textDecrypted.substr(_Pos, 1);

                if(monogramsText.find(currentBigram) == monogramsText.end()) 
                    monogramsText[currentBigram] = 1;
                else
                    monogramsText[currentBigram]++;
            }
            for (auto &_Elem : monogramsText)
                _Elem.second /= countAllSymbol;
     
            bigramsText.clear();
            for (int _Pos = 0; _Pos < textDecrypted.size() - 1; _Pos++) {
                std::wstring currentBigram = textDecrypted.substr(_Pos, 2);

                if(bigramsText.find(currentBigram) == bigramsText.end()) 
                    bigramsText[currentBigram] = 1;
                else
                    bigramsText[currentBigram]++;
            }
            for (auto &_Elem : bigramsText)
                _Elem.second /= countBigrams;

            trigramsText.clear();
            for (int _Pos = 0; _Pos < textDecrypted.size() - 2; _Pos++) {
                std::wstring currentTrigram = textDecrypted.substr(_Pos, 3);

                if(trigramsText.find(currentTrigram) == trigramsText.end()) 
                    trigramsText[currentTrigram] = 1;
                else
                    trigramsText[currentTrigram]++;
            }
            for (auto &_Elem : trigramsText)
                _Elem.second /= countTrigrams;


            _Efficiency = countDifference(monogramsText, monogramsFileMap, bigramsText, bigramsFileMap, trigramsText, trigramsFileMap);
            tempPopulationEfficiency[i].second = _Efficiency;
        }

        populationEfficiency = tempPopulationEfficiency;
    }


    std::sort(populationEfficiency.begin(), populationEfficiency.end(), 
        [] (const std::pair<std::wstring, double> &_Left, 
            const std::pair<std::wstring, double> &_Right) { return _Left.second < _Right.second; });


    monogramsKey = populationEfficiency[0].first;

    std::wstring bestKey;
    bestKey.resize(COUNT_RUS_LETTER);

    size_t _Pos = 0;
    for (const auto &x : monogramsKey) {
        bestKey[getIndexArray(fileKey[_Pos])] = x;
        _Pos++;
    }
    
    std::wcout << std::endl  << std::endl  << std::endl;
    std::wcout << bestKey << std::endl;

}


double metric(double _F, double _T, int64_t& _Count) {
    _Count++;
    return fabs(_F - _T) / std::max(_F, _T);
}

double Crypter::countDifference(std::unordered_map<std::wstring, double> &monogramsText,
                                std::unordered_map<std::wstring, double> &monogramsFile,
                                std::unordered_map<std::wstring, double> &bigramsText,
                                std::unordered_map<std::wstring, double> &bigramsFile,
                                std::unordered_map<std::wstring, double> &trigramsText,
                                std::unordered_map<std::wstring, double> &trigramsFile) {

    int64_t count1 = 0, count2 = 0, count3 = 0;;
    double result1 = 0, result2 = 0, result3 = 0;

    for (auto &x : monogramsText) {
        std::unordered_map<std::wstring, double>::const_iterator got = monogramsFile.find(x.first);
        if(got == monogramsFile.end()) {
            count1++;
            result1++;
        } else 
            result1 += metric(got->second, x.second, count1);
    }

    for (auto &x : bigramsText) {
        std::unordered_map<std::wstring, double>::const_iterator got = bigramsFile.find(x.first);
        if(got == bigramsFile.end()) {
            count2++;
            result2++;
        } else 
            result2 += metric(got->second, x.second, count2);
    }

    for (auto &x : trigramsText) {
        std::unordered_map<std::wstring, double>::const_iterator got = trigramsFile.find(x.first);
        if(got == trigramsFile.end()) {
            count3++;
            result3++;
        } else 
            result3 += metric(got->second, x.second, count3);
    }

    result1 /= count1;
    result2 /= count2;
    result3 /= count3;

    double weight1 = 1, weight2 = 2, weight3 = 3;
	double weightSum = weight1 + weight2 + weight3;

	return ( 
        weight1 * result1 + 
        weight2 * result2 + 
        weight3 * result3 
    ) / weightSum;
}



_CRYPTER_END

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

        void setTimeEnd();
        float getElapsed() const;

        friend std::ostream& operator<< (std::ostream& out, Timer& point) {
            return out << "SECONDS: " << point.getElapsed() << std::endl;
        }
    };
    void  Timer::setTimeEnd()       { timeEnd = high_resolution_clock::now(); }
    float Timer::getElapsed() const { return duration<float>(timeEnd - timeStart).count(); }
}

#undef _STD
#endif // _CRYPTER_H_