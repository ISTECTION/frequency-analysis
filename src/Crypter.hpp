#ifndef _CRYPTER_H_
#define _CRYPTER_H_
#include <iostream>
#include <string>
#include <codecvt>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <random>
#include <regex>
#include <locale>
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

const uint16_t COUNT_RUSSIAN_ABC = 33;

//? const PATH ------------------------------------------------------ ///
const _FILESYSTEM path ASSAY           = "texts/assay.txt";
const _FILESYSTEM path ENCRYPTED_ASSAY = "texts/encrypted_assay.txt";

const _FILESYSTEM path COMPOSITION = "texts/composition.txt";

const _FILESYSTEM path KEY = "output/key.txt";
const _FILESYSTEM path RUSSIAN_MONOGRAMS = "rus/russian_monograms.txt";
//? const PATH ------------------------------------------------------ ///

class Crypter
{
private:
    _STD wstring _Key;

    void locale() const;
    uint checkLetter(const uint) const;
    uint getLowerID_RU_UTF8(const uint) const;

    void openFileUTF8(_STD wstring &_Key, const _FILESYSTEM path); //? delete
    _STD wstring formatingText(_STD wstring &_Text);
    void codeConvertUTF8(std::wfstream &_Stream);

    uint getSeed() const;
    uint getRandom(std::mt19937 &engine, uint, uint) const;
public:
    Crypter() { 
        locale();
        // std::wstring key = generateRandomKey();

        std::wstring key = inKeyFromFileUTF8(KEY);
    
        std::wstring text = DecryptText(key, ENCRYPTED_ASSAY);
        wcout_t(text);

        // outFileUTF8(key, KEY);
        // outFileUTF8(text, ENCRYPTED_ASSAY);
    };
    ~Crypter() { };

    _STD wstring EncryptText(std::wstring& key, const _FILESYSTEM path);
    _STD wstring DecryptText(std::wstring& key, const _FILESYSTEM path);

    _STD wstring generateRandomKey() const;
    _STD wstring getRussianABC()     const;

    _STD wstring inKeyFromFileUTF8 (const _FILESYSTEM path);
    void outFileUTF8 (_STD wstring &_Text, const _FILESYSTEM path);
};

_STD wstring Crypter::EncryptText(std::wstring& key, const _FILESYSTEM path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) std::cerr << "File was not open";

    codeConvertUTF8(fin);
    std::wstring text, line;
    while (std::getline(fin, line)) 
        text += line;
    fin.close();

    std::wstring formattedText = formatingText(text);
    std::wstring ABC = getRussianABC();
    std::map<wchar_t, wchar_t> dict;

    for_t (ABC.size()) dict[ABC[_Pos]] = key[_Pos];
    for_t (formattedText.size()) 
        formattedText[_Pos] = dict[formattedText[_Pos]];

    return formattedText;
}

_STD wstring Crypter::DecryptText(std::wstring& key, const _FILESYSTEM path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    if(!fin) std::cerr << "File was not open";
    codeConvertUTF8(fin);
    std::wstring text, line;
    while (std::getline(fin, line)) 
        text += line;
    fin.close();

    std::wstring ABC = getRussianABC();
    std::map<wchar_t, wchar_t> dict;

    for_t (ABC.size()) dict[key[_Pos]] = ABC[_Pos];
    for_t (text.size()) 
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

void Crypter::openFileUTF8(_STD wstring& _Key, const _FILESYSTEM path _Path) {
    std::wfstream fin(_Path, std::wfstream::in); 
    codeConvertUTF8(fin);
    for (std::wstring line; std::getline(fin, line);) {
        std::wstringstream wstream(line);
        wchar_t _Symbol;
        wstream >> _Symbol;
        _Key.push_back(getLowerID_RU_UTF8(static_cast<uint>(_Symbol)));
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

_CRYPTER_END
#undef _STD
#endif // _CRYPTER_H_