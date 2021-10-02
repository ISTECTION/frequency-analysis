#include <iostream>
#include <string>
#include <filesystem>

#include "Crypter.hpp"
using namespace crypter;

int main(int argc, char *argv[])
{
    Crypter A;
    // std::wstring key = A.getRussianABC();

    // std::wstring text = A.EncryptText(key, ); 



    // A.outFileUTF8(text, std::filesystem::path("texts/Jules_Verne_encrypted.txt"));


    // std::wstring key = A.getRussianABC();
    // std::wcout << key;
    // std::wstring text = A.EncryptText(key, COMPOSITION);


    // A.outFileUTF8(key, KEY);
    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);

    // A.FrequencyAnalysis(COMPOSITION_ENCRYPTED);

    // std::wstring key = A.generateRandomKey();
    // A.outFileUTF8(key, KEY);
    // std::wstring text = A.EncryptText(key, COMPOSITION);
    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);


    std::wstring key = A.inKeyFromFileUTF8(KEY);
    std::wstring text = A.DecryptText(key, COMPOSITION_ENCRYPTED);

    std::wcout << text;


    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);

    // wcout_t(text);

    // A.outFileUTF8(key, KEY);
    // outFileUTF8(text, ENCRYPTED_ASSAY);
}