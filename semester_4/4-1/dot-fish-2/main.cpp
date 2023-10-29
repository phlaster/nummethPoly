#include "./headers/routines.hpp"

int main(int argc, char** argv){  
    if (argc == 1){
        cout << "flags to execute tasks:\n"
            << "\t--all\n"
            << "\t--baza\n"
            << "\t--minimum / --min\n"
            << "\t--dostatochno / --dost\n";
        return 0;
    }
    bool 
        do_baza = false,
        do_minimum = false,
        do_dostatochno = false;

    for(int i = 1; i < argc; i++){
        string s = argv[i];
        if (!s.compare("--all")){
            Baza();
            Minimum();
            Dostatochno();
            return 0;
        }
        else if (!s.compare("--baza")){
            do_baza = true;
        }
        else if (!(s.compare("--min")) || !(s.compare("--minimum"))){
            do_minimum = true;
        }
        else if (!(s.compare("--dostatochno")) || !(s.compare("--dost"))){
            do_dostatochno = true;
        }
    }
    if (do_baza)
        Baza();
    if (do_minimum)
        Minimum();
    if (do_dostatochno)
        Dostatochno();

    return 0;
}