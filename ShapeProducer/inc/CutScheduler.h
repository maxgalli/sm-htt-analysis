#ifndef CUTSCHEDULER
#define CUTSCHEDULER

using namespace ROOT;

class CutScheduler {

    public:
        // Default constructor
        CutScheduler();

        // Only name passed
        CutScheduler(std::string name);

        // Name and vector string
        CutScheduler(std::string name, std::vector<Cut> cuts);

        CutScheduler(std::string name, RDataFrame rdf, std::vector<Cut> cuts);

        stringvec getCuts(std::string name);
        void setCuts(stringvec v);
        void applyCuts();
        RDataFrame getRDataFrame() {return m_rdf;}

    private:
        RDataFrame m_rdf = RDataFrame(0);
        std::string m_name;
        stringvec m_cut_names_vec;
        std::vector<Cut> m_cuts_vec;
};

ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> ScheduleCuts(RDataFrame &, std::vector<Cut>);

#endif

