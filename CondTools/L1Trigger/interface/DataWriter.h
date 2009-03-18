#ifndef CondTools_L1Trigger_DataWriter_h
#define CondTools_L1Trigger_DataWriter_h

// Framework
#include "FWCore/Framework/interface/IOVSyncValue.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/DataKey.h"

#include "CondCore/DBCommon/interface/TypedRef.h"
#include "CondCore/MetaDataService/interface/MetaData.h"

#include "DataFormats/Provenance/interface/RunID.h"

// L1T includes
#include "CondFormats/L1TObjects/interface/L1TriggerKeyList.h"
#include "CondFormats/L1TObjects/interface/L1TriggerKey.h"

#include "CondTools/L1Trigger/interface/WriterProxy.h"

#include <string>
#include <map>

namespace l1t
{

/* This class is used to write L1 Trigger configuration data to Pool DB.
 * It also has a function for reading L1TriggerKey directly from Pool.
 *
 * In order to use this class to write payloads, user has to make sure to register datatypes that she or he is
 * interested to write to the framework. This should be done with macro REGISTER_L1_WRITER(record, type) found in
 * WriterProxy.h file. Also, one should take care to register these data types to CondDB framework with macro
 * REGISTER_PLUGIN(record, type) from registration_macros.h found in PluginSystem.
 */

class DataWriter
{
 public:
  explicit DataWriter() {} ;
  virtual ~DataWriter () {};

  // Payload and IOV writing functions.  

  // Get payload from EventSetup and write to DB with no IOV
  // recordType = "record@type", return value is payload token
  std::string writePayload( const edm::EventSetup& setup,
			    const std::string& recordType ) ;

  // Use PoolDBOutputService to append IOV with sinceRun to IOV sequence
  // for given ESRecord.  PoolDBOutputService knows the corresponding IOV tag.
  // Return value is true if IOV was updated; false if IOV was already
  // up to date.
  bool updateIOV( const std::string& esRecordName,
		  const std::string& payloadToken,
		  edm::RunNumber_t sinceRun ) ;

  // Write L1TriggerKeyList payload and set IOV.  Takes ownership of pointer.
  void writeKeyList( L1TriggerKeyList* keyList,
		     edm::RunNumber_t sinceRun = 0 ) ;

  // Read L1TriggerKey directly from Pool, not from EventSetup.
  void readKey( const std::string& payloadToken,
		L1TriggerKey& outputKey ) ;


 protected:
};

} // ns

#endif
