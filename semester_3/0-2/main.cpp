// Вариант 11 Б
/*
    f1(x) = ctg(x) + x^2
    f2(x) = x^5 - 3.2x^3 + 2.5x^2 - 7x + 1.5
*/
#include "functions.hpp"
#include <thread>

const Vec LIMS_1 = {0.5, 2.75};
const Vec LIMS_2 = {-2.4, 2.1};


void do_task(int numer, double (*f)(double, bool), int minNodes, int maxNodes, int coeffMult, double distModulo=0)
{
    Vec lims;
    Str threadname;
    Buffer buffer;
    Str summaryName = "CSVs/task";

    if (f==f1){
        lims = LIMS_1;
        threadname = "f1";
        
    } else {
        lims = LIMS_2;
        threadname = "f2";
    }
    switch (numer)
    {
        case 3:
        {
            task2_3(f, minNodes, maxNodes, coeffMult, lims, threadname, ref(buffer));
            summaryName +="2-3_"+threadname+"_summ_"+to_string(minNodes)+"-"+to_string(maxNodes);
            break;
        }

        case 4:
        {
            task4(f, minNodes, maxNodes, coeffMult, lims, threadname, ref(buffer));
            summaryName += "4_"+threadname+"_summ_"+to_string(minNodes)+"-"+to_string(maxNodes);
            break;
        }

        case 5:
        {
            task5(f, minNodes, maxNodes, coeffMult, distModulo, lims, threadname, ref(buffer));
            summaryName += "5_"+threadname+"_summ_"+to_string(minNodes)+"-"+to_string(maxNodes)+"_"+to_string(distModulo);
            break; 
        }
        default:
            cerr << "Wrong Task number " << numer << endl;
        break;
    }

    ofstream outSummary(summaryName+"_.csv");
    if (!outSummary.is_open())
        throw runtime_error("Не удалось открыть файл для записи сводки!\n");
    outSummary << buffer.data;
    outSummary.close();
}

int main(int argc, char *argv[])
{
    int minNodes = stoi(argv[1]);
    int maxNodes = stoi(argv[2]);
    int interpolCoeff = stoi(argv[3]);
    double distModulo = fabs(stod(argv[4]));
    assert(minNodes <= maxNodes && "Сначала нижняя граница!\n");
    thread t1, t2;

    for (int i=3; i<=5; i++)
    {
        t1 = thread(do_task, i, f1, minNodes, maxNodes, interpolCoeff, distModulo);
        t2 = thread(do_task, i, f2, minNodes, maxNodes, interpolCoeff, distModulo);
        t1.join(); t2.join();
        cout <<  "Задание " << i << " выполнено!\n";
    }
    return 0;
}
