#include "bubble/model.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <stdexcept>
int main(int argc,char**argv){if(argc!=5)throw std::runtime_error("usage: input checkpoint config output next_dt");
  std::ifstream in(argv[1]);auto root=nlohmann::json::parse(in);auto c=bubble::Config::load(argv[2]);auto p=bubble::make_checkpoint_provenance(c);
  if(root.at("provenance").at("pore_sha256")!=p.pore_sha256||root.at("provenance").at("throat_sha256")!=p.throat_sha256||root.at("provenance").at("connect_sha256")!=p.connect_sha256)
    throw std::runtime_error("network hashes differ; retag prohibited");
  root["provenance"]["continuation_config_sha256"]=p.continuation_config_sha256;root["next_dt_s"]=std::stod(argv[4]);root["consecutive_fast_steps"]=0;
  std::ofstream out(argv[3]);out<<std::setprecision(17)<<root.dump(1)<<'\n';}
