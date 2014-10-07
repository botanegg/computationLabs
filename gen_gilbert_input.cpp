#include <iostream>
#include <fstream>

using namespace std;

int main() {
    fstream file("lab02/input.txt");
    int size = 0;
    cin >> size;
    int size_d1 = size - 1;

    file << size << endl;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file << 1.0 / (i + j + 1);

            if (j < size_d1) file << ' ';
        }
        file << endl;
    }

    file.close();

    return 0;
}