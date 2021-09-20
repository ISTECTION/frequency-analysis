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
#define for_t(_Size) for(size_t _Pos = 0; _Pos < (_Size); _Pos++)
#define wcout_t(_Line) std::wcout << _Line << std::endl

#define _STD        ::std::
#define _FILESYSTEM ::std::filesystem::

_CRYPTER_BEGIN
using uint = unsigned int;
using SCS = std::chrono::system_clock;
using pair_wchar = std::pair<wchar_t, uint>;


const uint16_t COUNT_RUSSIAN_ABC = 33;

//? const PATH ------------------------------------------------------ ///
const _FILESYSTEM path KEY = "output/key.txt";

const _FILESYSTEM path ASSAY = "texts/assay.txt";
const _FILESYSTEM path ASSAY_ENCRYPTED = "texts/assay_encrypted.txt";

const _FILESYSTEM path COMPOSITION = "texts/composition.txt";
const _FILESYSTEM path COMPOSITION_ENCRYPTED = "texts/composition_encrypted.txt";

const _FILESYSTEM path RUSSIAN_MONOGRAMS = "rus/russian_monograms.txt";
const _FILESYSTEM path RUSSIAN_BIGRAMS = "rus/russian_bigrams.txt";
//? const PATH ------------------------------------------------------ ///

/// Index: 6 - Russian `E
uint getIndexArray(const uint _Off) { return _Off < 1078 ? _Off - 1072 : _Off == 1105 ? 6 : _Off - 1071; }

class Crypter
{
private:
    void locale() const;
    uint checkLetter(const uint) const;
    _STD wstring getLowerStr(_STD wstring) const;
    uint getLowerID_RU_UTF8(const uint) const;

    _STD wstring formatingText(_STD wstring &_Text);
    void codeConvertUTF8(std::wfstream &_Stream);

    uint getSeed() const;
    uint getRandom(std::mt19937 &engine, uint, uint) const;

    void inFileUTF8_FORMAT_WCHAR_COUNT(std::vector<std::pair<std::wstring, uint>> &_V, const _FILESYSTEM path);
    double dispersion(std::unordered_map<std::wstring, uint> &_T, std::unordered_map<std::wstring, uint> &_F);
public:
    Crypter() { locale(); };
    ~Crypter() { };

    _STD wstring EncryptText(std::wstring& key, const _FILESYSTEM path);
    _STD wstring DecryptText(std::wstring& key, const _FILESYSTEM path);

    _STD wstring generateRandomKey() const;
    _STD wstring getRussianABC()     const;

    _STD wstring inKeyFromFileUTF8 (const _FILESYSTEM path);
    void outFileUTF8 (_STD wstring &_Text, const _FILESYSTEM path);

    void FrequencyAnalysis(const _FILESYSTEM path);
};

_STD wstring Crypter::EncryptText(std::wstring& key, const _FILESYSTEM path _Path) {
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

    for_t (ABC.size()) dict[ABC[_Pos]] = key[_Pos];

    for(size_t _Pos = 0; _Pos < formattedText.size(); _Pos++)
        formattedText[_Pos] = dict[formattedText[_Pos]];

    return formattedText;
}

_STD wstring Crypter::DecryptText(std::wstring& key, const _FILESYSTEM path _Path) {
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

    for_t (ABC.size()) dict[key[_Pos]] = ABC[_Pos];

    for(size_t _Pos = 0; _Pos < text.size(); _Pos++)
        text[_Pos] = dict[text[_Pos]];

    return text;
}

_STD wstring Crypter::formatingText(_STD wstring &_Text) {
    for (auto &pos : _Text) 
        pos = checkLetter(pos);

    std::wcmatch charRg;    
    std::wregex rg (L"[а-яё]");
    std::wstring formattedText;
    while (std::regex_search(_Text.c_str(), charRg, rg)) {
        formattedText += charRg.str();
        _Text = charRg.suffix().str();
    }
    return formattedText;
}

_STD wstring Crypter::generateRandomKey() const {
    std::mt19937 engine;
    engine.seed(getSeed());

    std::wstring ABC = getRussianABC();
    std::wstring key;
    key.resize(COUNT_RUSSIAN_ABC);

    for(size_t _Pos = 0; _Pos < COUNT_RUSSIAN_ABC; _Pos++) {
        uint rdPos = getRandom(engine, ZERO, COUNT_RUSSIAN_ABC - _Pos);
        key[_Pos] = ABC[rdPos];
        ABC.erase(rdPos, 1);
    }
    return key;
}

_STD wstring Crypter::getRussianABC() const {
    const size_t LETTER_A = 1072, INC_LETTER_E = 1078;
    std::wstring ABC;
    ABC.resize(COUNT_RUSSIAN_ABC);

    for (size_t _Pos = 0; _Pos < INC_LETTER_E - LETTER_A; _Pos++) 
        ABC[_Pos] = LETTER_A + _Pos;
    ABC[INC_LETTER_E - LETTER_A] = 1105;

    for (size_t _Pos = INC_LETTER_E - LETTER_A; _Pos < COUNT_RUSSIAN_ABC - 1; _Pos++) 
        ABC[_Pos + 1] = LETTER_A + _Pos;

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
    std::wstring lowerStr;
    for(const auto &_Symbol : _Str)
        lowerStr += getLowerID_RU_UTF8(_Symbol);
    
    return lowerStr;
}

void Crypter::inFileUTF8_FORMAT_WCHAR_COUNT(std::vector<std::pair<std::wstring, uint>> &_V,
     const _FILESYSTEM path _Path) {

    std::wfstream fin(_Path, std::wfstream::in);
    if (!fin) 
        std::cerr << "File was not open";
    codeConvertUTF8(fin);
    std::wstring line;
    while (std::getline(fin, line)) {
        std::wstringstream wstream(line);
        std::wstring _Str;
        uint _Size;
        wstream >> _Str >> _Size;
        _V.push_back(std::pair(getLowerStr(_Str), _Size));
    }
    fin.close();
}

_STD wstring Crypter::inKeyFromFileUTF8(const _FILESYSTEM path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if (!fin) std::cerr << "File was not open";
    codeConvertUTF8(fin);
    std::wstring key;
    std::getline(fin, key);
    fin.close();
    return key;
}

void Crypter::outFileUTF8(_STD wstring& _Text, const _FILESYSTEM path _Path) {
    std::wfstream fout(_Path, std::wfstream::out);
    codeConvertUTF8(fout);
    fout << _Text;
    fout.close();
}

void Crypter::codeConvertUTF8(std::wfstream &_Stream) { 
    _Stream.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>)); 
}

void Crypter::FrequencyAnalysis(const _FILESYSTEM path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) 
        std::cerr << "File was not open";

    codeConvertUTF8(fin);
    std::wstring text, line;
    while (std::getline(fin, line)) 
        text += line;
    fin.close();

    std::vector<std::pair<wchar_t, uint>> countMonograms;
    countMonograms.resize(COUNT_RUSSIAN_ABC);

    std::wstring ABC = getRussianABC();
    for(size_t _Pos = 0; _Pos < ABC.size(); _Pos++) {
        countMonograms[_Pos].first = ABC[_Pos];
        countMonograms[_Pos].second = ZERO;
    }

    for (const auto &_El : text) 
        countMonograms[getIndexArray(_El)].second++;
    
    uint countAllSymbol = std::accumulate(countMonograms.begin(), countMonograms.end(), 0, 
        [](const std::size_t _S, const auto &_Elem) { return _S + _Elem.second; });

    std::sort(countMonograms.begin(), countMonograms.end(), 
        [] (const pair_wchar &_Left, const pair_wchar &_Right) { return _Left.second > _Right.second; });

    std::wstring monogramsKey;
    std::wostringstream w_stream;
    for(auto &_El : countMonograms)
        w_stream << _El.first;
    monogramsKey = w_stream.str();

    std::vector<std::pair<std::wstring, uint>> monogramsFile;
    inFileUTF8_FORMAT_WCHAR_COUNT(monogramsFile, RUSSIAN_MONOGRAMS);

    std::unordered_map<wchar_t, wchar_t> monogramsMap;
    for (size_t _Pos = 0; _Pos < monogramsKey.size(); _Pos++)
        monogramsMap[monogramsKey[_Pos]] = monogramsFile[_Pos].first.at(0);

    std::wstring textDecrypted;
    for (auto &x : text) 
        textDecrypted += monogramsMap[x];

    std::unordered_map<std::wstring, uint> bigramsText;
    for (int _Pos = 0; _Pos < textDecrypted.size() - 1; _Pos++) {
        std::wstring currentBigram = textDecrypted.substr(_Pos, 2);

        if(bigramsText.find(currentBigram) == bigramsText.end()) 
            bigramsText[currentBigram] = 1;
        else
            bigramsText[currentBigram]++;
    }
    std::vector<std::pair<std::wstring, uint>> bigramsFile;
    inFileUTF8_FORMAT_WCHAR_COUNT(bigramsFile, RUSSIAN_BIGRAMS);
    std::unordered_map<std::wstring, uint> bigramsFileMap;
    for (const auto &x : bigramsFile)
        bigramsFileMap[x.first] = x.second;

    std::wstring reservMonogramsKey = monogramsKey;
    double _Point = dispersion(bigramsText, bigramsFileMap);
    std::cout << _Point << std::endl;


    for (size_t _Off = 0; _Off < COUNT_RUSSIAN_ABC - 1; _Off++) {
        std::swap(monogramsKey[_Off], monogramsKey[_Off + 1]);

        monogramsMap.clear();
        for (size_t _Pos = 0; _Pos < monogramsKey.size(); _Pos++)
            monogramsMap[monogramsKey[_Pos]] = monogramsFile[_Pos].first.at(0);

        textDecrypted.clear();
        for (auto &x : text) 
            textDecrypted += monogramsMap[x];

        bigramsText.clear();
        for (int _Pos = 0; _Pos < textDecrypted.size() - 1; _Pos++) {
            std::wstring currentBigram = textDecrypted.substr(_Pos, 2);

            if(bigramsText.find(currentBigram) == bigramsText.end()) 
                bigramsText[currentBigram] = 1;
            else
                bigramsText[currentBigram]++;
        }

        double current_Point = dispersion(bigramsText, bigramsFileMap);
        std::cout << current_Point << std::endl;

        if (_Point < current_Point)
            std::swap(monogramsKey[_Off], monogramsKey[_Off + 1]);
        else _Point = current_Point;
    }

    std::wcout << textDecrypted;

    // for(const auto &x : bigramsFile)
    //     std::wcout << x.first << " " << x.second << std::endl;
}

double Crypter::dispersion(std::unordered_map<std::wstring, uint> &_T, std::unordered_map<std::wstring, uint> &_F) {
    double _Point = 0;
    for (const auto &_Item : _T) 
        if (_F.find(_Item.first) != _F.end())
            _Point += pow(static_cast<double>(_Item.second) / _F[_Item.first], 2);     
    
    return _Point;
}

_CRYPTER_END
#undef _STD
#endif // _CRYPTER_H_