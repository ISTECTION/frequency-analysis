#include <iostream>
#include <string>
#include <filesystem>

#include "Crypter.hpp"
using namespace crypter;

int main(int argc, char *argv[])
{
    Crypter A;
    // std::wstring key = A.generateRandomKey();

    // std::wstring text = A.EncryptText(key, std::filesystem::path("texts/Jules_Verne.txt"));
    // A.outFileUTF8(text, std::filesystem::path("texts/Jules_Verne_encrypted.txt"));

    A.FrequencyAnalysis(std::filesystem::path("texts/Jules_Verne_encrypted.txt"));

    // std::wstring key = A.generateRandomKey();

    // std::wstring key = inKeyFromFileUTF8(KEY);

    // std::wstring text = A.EncryptText(key, COMPOSITION);
    // A.outFileUTF8(text, COMPOSITION_ENCRYPTED);

    // wcout_t(text);

    // A.outFileUTF8(key, KEY);
    // outFileUTF8(text, ENCRYPTED_ASSAY);
}