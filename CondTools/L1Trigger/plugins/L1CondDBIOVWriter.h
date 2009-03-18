#ifndef CondTools_L1Trigger_L1CondDBIOVWriter_h
#define CondTools_L1Trigger_L1CondDBIOVWriter_h
// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1CondDBIOVWriter
// 
/**\class L1CondDBIOVWriter L1CondDBIOVWriter.h CondTools/L1Trigger/interface/L1CondDBIOVWriter.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Werner Sun
//         Created:  Sun Mar  2 20:10:36 CET 2008
// $Id: L1CondDBIOVWriter.h,v 1.2 2008/09/27 02:38:19 wsun Exp $
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondTools/L1Trigger/interface/DataWriter.h"

// forward declarations

class L1CondDBIOVWriter : public edm::EDAnalyzer {
   public:
      explicit L1CondDBIOVWriter(const edm::ParameterSet&);
      ~L1CondDBIOVWriter();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      l1t::DataWriter m_writer ;
      std::string m_tscKey ;

      // List of record@type, used only for objects not tied to TSC key.
      // Otherwise, list of records comes from L1TriggerKey.
      std::vector< std::string > m_recordTypes ;

      // When true, set IOVs for objects not tied to the TSC key.  The records
      // and objects to be updated are given in the toPut parameter, and
      // m_tscKey is taken to be a common key for all the toPut objects, not
      // the TSC key.  The IOV for L1TriggerKey is not updated when
      // m_ignoreTriggerKey = true.
      bool m_ignoreTriggerKey ;
};


#endif
