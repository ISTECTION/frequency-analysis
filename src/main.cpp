#include <iostream>
#include <string>
#include <filesystem>

#include "Crypter.hpp"
#include "Utils.hpp"


using namespace crypter;
using namespace utils;

int main(int argc, char *argv[])
{
///TODO:
    Crypter A;
///TODO:

    // std::wstring key = A.generateRandomKey();
    // std::wstring text = A.EncryptText(key, "texts/gogol.txt");

    // A.outFileUTF8(key, "output/gogol_key.txt");
    // A.outFileUTF8(text, "texts/gogol_encrypted.txt");


    // Utils a;
    // a.MONO();
    // a.BI();
    // a.TRI();

    // std::wstring key = A.getRussianABC();
    // std::wstring text = A.EncryptText(key, ); 


    // A.outFileUTF8(text, std::filesystem::path("texts/Jules_Verne_encrypted.txt"));
    // std::wstring key = A.getRussianABC();
    // std::wcout << key;
    // std::wstring text = A.EncryptText(key, COMPOSITION);
    // A.outFileUTF8(key, KEY);
    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);


///TODO:
    Timer::Timer t;
    std::wstring key = A.FrequencyAnalysis(COMPOSITION_ENCRYPTED);
    t.setTimeEnd();
    std::cout <<  std::endl <<  t << std::endl;
    std::wcout << std::endl << key;
///TODO:



    // std::wstring key = A.generateRandomKey();
    // A.outFileUTF8(key, KEY);
    // std::wstring text = A.EncryptText(key, COMPOSITION);
    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);
    // std::wstring key = A.inKeyFromFileUTF8(KEY);


///TODO:
    std::wstring text = A.DecryptText(key, COMPOSITION_ENCRYPTED);
    std::wcout << std::endl << std::endl << text;
///TODO:


    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);
    // wcout_t(text);
    // A.outFileUTF8(key, KEY);
    // outFileUTF8(text, ENCRYPTED_ASSAY);
}