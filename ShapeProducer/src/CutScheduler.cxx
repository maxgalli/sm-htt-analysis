#include "CutScheduler.h"

CutScheduler::CutScheduler(std::string name, RDataFrame rdf, std::vector<Cut> cuts_vec) {
    m_name = name;
    m_cuts_vec = cuts_vec;
    m_rdf = rdf;
}

void CutScheduler::applyCuts() {
    for (auto cut : m_cuts_vec) m_rdf.Filter(cut.getStringcut());
}
