#include "CutScheduler.h"

CutScheduler::CutScheduler(std::string name, RDataFrame rdf, std::vector<Cut> cuts_vec) {
    m_name = name;
    m_cuts_vec = cuts_vec;
    m_rdf = rdf;
}

void CutScheduler::applyCuts() {
    for (auto cut : m_cuts_vec) m_rdf.Filter(cut.getStringcut());
}

ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> ScheduleCuts(RDataFrame& rdf, std::vector<Cut> cuts_vec) {
    int len = cuts_vec.size();
    std::vector<ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>> rdf_vec;
    rdf_vec.push_back(
            rdf.Filter(
                cuts_vec[0].getStringcut()
                )
            );
    int rdfvec_counter = 0;
    for (auto cut : cuts_vec) {
        rdf_vec.push_back(
                rdf_vec[rdfvec_counter].
                Filter(cut.getStringcut()
                    )
                );
        rdfvec_counter++;
    }
    return rdf_vec.back();
}
