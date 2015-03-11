// Ver 1.2

#include <string>
#include <vector>

#ifndef EOS_h
#define EOS_h


class EOS
{
    private:
        std::vector<double>* ed_table;
        std::vector<double>* p_table;
        std::vector<double>* sd_table;
        std::vector<double>* T_table;
        long table_length;
        double delta_ed, max_ed;
        double p1,p2,s1,s2,T1,T2;
    public:
        EOS();
        EOS(char*);
        EOS(char*,char*);
        void loadEOSFromFile(char*);
        void loadEOSFromFile(char*,char*);
        double p(double ed);
        double sd(double ed);
        double T(double ed);
        double edFromP(double sd);
        double edFromSd(double sd);
        double edFromT(double sd);
};

#endif

// Ver 1.2:
// -- Varible max_ed variable added.
// -- Variable p1,p2,s1,s2,T1,T2 added (for EOS extrapolation).
// -- The EOS and loadEOSFromFile functions are overloaded to provide
//    possible support for extrapolation using coeff file.
