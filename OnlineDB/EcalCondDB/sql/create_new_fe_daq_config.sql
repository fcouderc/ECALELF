/* RUN TYPE DEFINITIONS */

 
/* THIS IS THE MASTER COND2CONF TABLE */ 

CREATE TABLE COND2CONF_TYPE_DEF (
	DEF_ID NUMBER(2) NOT NULL,
	REC_TYPE VARCHAR2(20)
);

ALTER TABLE cond2conf_type_def ADD CONSTRAINT COND2CONF_TYPE_DEF_PK PRIMARY KEY (DEF_ID);

INSERT INTO COND2CONF_TYPE_DEF(DEF_ID, REC_TYPE) VALUES( '0', 'PEDESTAL_OFFSETS' );
INSERT INTO COND2CONF_TYPE_DEF(DEF_ID, REC_TYPE) VALUES( '1', 'DELAY_OFFSETS'    );
INSERT INTO COND2CONF_TYPE_DEF(DEF_ID, REC_TYPE) VALUES( '2', 'DCC_WEIGHTS'  );
INSERT INTO COND2CONF_TYPE_DEF(DEF_ID, REC_TYPE) VALUES( '3', 'BAD_CRYSTALS'  );
INSERT INTO COND2CONF_TYPE_DEF(DEF_ID, REC_TYPE) VALUES( '4', 'BAD_TT'  );

CREATE TABLE COND2CONF_INFO (
  REC_ID		NUMBER(10) NOT NULL,
  REC_TYPE_ID           NUMBER(10) ,          
  REC_date              DATE ,
  LOCATION_ID           NUMBER(10),
  RUN_NUMBER  NUMBER(10),
  short_desc         VARCHAR2(100),
  db_timestamp          TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE cond2conf_info ADD CONSTRAINT COND2CONF_INFO_pk PRIMARY KEY (REC_ID);
ALTER TABLE COND2CONF_INFO ADD CONSTRAINT COND2CONF_INFO_FK FOREIGN KEY (REC_TYPE_ID) REFERENCES COND2CONF_TYPE_DEF(DEF_ID);

CREATE SEQUENCE COND2CONF_INFO_SQ INCREMENT BY 1 START WITH 1;

/* PEDESTAL OFFSETS */

CREATE TABLE pedestal_offsets_info (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE pedestal_offsets_INFO ADD CONSTRAINT  pedestal_offsets_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE pedestal_offsets_INFO ADD CONSTRAINT pedestal_offsets_INFO_uk UNIQUE(tag,version);


CREATE TABLE pedestal_offsets_dat (
REC_id NUMBER(10) NOT NULL,
sm_id  NUMBER(10),
fed_id NUMBER(10),
tt_id  NUMBER(10),
cry_id NUMBER(10), 
low    NUMBER, 
mid    NUMBER, 
high   NUMBER
);

ALTER TABLE PEDESTAL_OFFSETS_DAT ADD CONSTRAINT PEDESTAL_OFFSETS_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE PEDESTAL_OFFSETS_DAT ADD CONSTRAINT PEDESTAL_OFFSETS_DAT_pk PRIMARY KEY (rec_id, sm_id,tt_id, cry_id );


/* THIS IS FOR THE DELAYS */ 


CREATE TABLE DELAYS_info (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE delays_INFO ADD CONSTRAINT  delays_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE delays_INFO ADD CONSTRAINT delays_INFO_uk UNIQUE(tag,version);


CREATE TABLE DELAYS_DAT (
  REC_ID		NUMBER(10) NOT NULL,
  SM_ID NUMBER(10),
  FED_ID NUMBER(10),
  TT_ID NUMBER(10),
  TIME_OFFSET NUMBER
);

ALTER TABLE DELAYS_DAT ADD CONSTRAINT DELAYS_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE DELAYS_DAT ADD CONSTRAINT DELAYS_DAT_pk PRIMARY KEY (rec_id, sm_id,tt_id);



/* THIS IS FOR THE WEIGHTS */ 
/* THIS IS FOR THE WEIGHTS */ 

CREATE TABLE dcc_weights_info (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);
ALTER TABLE dcc_weights_INFO ADD CONSTRAINT dcc_weights_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE dcc_weights_INFO ADD CONSTRAINT dcc_weights_INFO_uk UNIQUE(tag,version);


CREATE TABLE DCC_WEIGHTS_DAT (
  REC_ID                NUMBER(10) NOT NULL,
  SM_ID NUMBER(10),
  FED_ID NUMBER(10),
  TT_ID NUMBER(10),
  CRY_ID NUMBER(10),
  WEI0 NUMBER,
  WEI1 NUMBER,
  WEI2 NUMBER,
  WEI3 NUMBER,
  WEI4 NUMBER,
  WEI5 NUMBER
);
ALTER TABLE DCC_WEIGHTS_DAT ADD CONSTRAINT DCC_WEIGHTS_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE DCC_WEIGHTS_DAT ADD CONSTRAINT DCC_WEIGHTS_DAT_pk PRIMARY KEY (rec_id, sm_id,tt_id, cry_id );


CREATE TABLE DCC_WEIGHTSAMPLE_DAT (
  REC_ID NUMBER(10) NOT NULL,
  FED_ID NUMBER(10),
  sample_id NUMBER(2),
  weight_number NUMBER(1)
);
ALTER TABLE DCC_WEIGHTSAMPLE_DAT ADD CONSTRAINT DCC_WEIGHTSAMPLE_DAT_pk PRIMARY KEY (rec_id, fed_id, sample_id );
ALTER TABLE DCC_WEIGHTSAMPLE_DAT ADD CONSTRAINT DCC_WEIGHTSAMPLE_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_
ID);


create OR REPLACE TRIGGER cond2conf_auto_tg3
  BEFORE INSERT ON dcc_weights_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('DCC_WEIGHTS', :newiov.rec_id, :newiov.tag, :newiov.version)
/



/* CRYSTAL bad channels */

CREATE TABLE BAD_Crystals_info (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE bad_crystals_INFO ADD CONSTRAINT bad_crystals_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE bad_crystals_INFO ADD CONSTRAINT bad_crystals_INFO_uk UNIQUE(tag,version);


CREATE TABLE bad_crystals_dat (
REC_id NUMBER(10) NOT NULL,
sm_id  NUMBER(10),
fed_id NUMBER(10),
tt_id  NUMBER(10),
cry_id NUMBER(10), 
status    NUMBER 
);

ALTER TABLE bad_crystals_DAT ADD CONSTRAINT bad_crystals_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE bad_crystals_DAT ADD CONSTRAINT bad_crystals_DAT_pk PRIMARY KEY (rec_id, sm_id,tt_id, cry_id );

/* TT bad channels */

CREATE TABLE BAD_TT_info (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE bad_tt_INFO ADD CONSTRAINT bad_tt_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE bad_tt_INFO ADD CONSTRAINT bad_tt_INFO_uk UNIQUE(tag,version);


CREATE TABLE bad_tt_dat (
REC_id NUMBER(10) NOT NULL,
sm_id  NUMBER(10),
fed_id NUMBER(10),
tt_id  NUMBER(10),
status    NUMBER 
);

ALTER TABLE bad_tt_DAT ADD CONSTRAINT bad_tt_DAT_FK FOREIGN KEY (REC_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE bad_tt_DAT ADD CONSTRAINT bad_tt_DAT_pk PRIMARY KEY (rec_id, sm_id,tt_id );



/* towers to bypass */

CREATE TABLE TOWERS_TO_BYPASS_INFO (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE TOWERS_TO_BYPASS_INFO ADD CONSTRAINT TOWERS_TO_BYPASS_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE TOWERS_TO_BYPASS_INFO ADD CONSTRAINT TOWERS_TO_BYPASS_INFO_uk UNIQUE(tag,version);


CREATE TABLE TOWERS_TO_BYPASS_dat (
REC_id NUMBER(10) NOT NULL,
fed_id NUMBER(3) not null,
tr_id  NUMBER(3) not null,
tt_id  NUMBER(3) not null,
time_corr NUMBER(10), 
status    NUMBER(10) 
);

ALTER TABLE TOWERS_TO_BYPASS_DAT ADD CONSTRAINT TOWERS_TO_BYPASS_DAT_pk PRIMARY KEY (rec_id, fed_id, tr_id, tt_id );

/* vfes to reject */

CREATE TABLE  VFES_TO_REJECT_INFO (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE VFES_TO_REJECT_INFO ADD CONSTRAINT VFES_TO_REJECT_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE VFES_TO_REJECT_INFO ADD CONSTRAINT VFES_TO_REJECT_INFO_uk UNIQUE(tag,version);


CREATE TABLE VFES_TO_REJECT_dat (
REC_id NUMBER(10) NOT NULL,
fed_id NUMBER(3) not null,
tt_id  NUMBER(3) not null,
vfe_id  NUMBER(3) not null,
gain NUMBER(10),
status    NUMBER(10)
);

ALTER TABLE VFES_TO_REJECT_DAT ADD CONSTRAINT VFES_TO_REJECT_DAT_pk PRIMARY KEY (rec_id, fed_id, tt_id, vfe_id );


/*  GOL bias current */

CREATE TABLE  GOL_BIAS_CURRENT_INFO (
 rec_id NUMBER(10) NOT NULL,
 TAG VARCHAR2(100),
 version NUMBER(10),
 db_timestamp  TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE GOL_BIAS_CURRENT_INFO ADD CONSTRAINT GOL_BIAS_CURRENT_INFO_PK PRIMARY KEY (rec_id);
ALTER TABLE GOL_BIAS_CURRENT_INFO ADD CONSTRAINT GOL_BIAS_CURRENT_INFO_uk UNIQUE(tag,version);

CREATE TABLE GOL_BIAS_CURRENT_DAT (
REC_id NUMBER(10) NOT NULL,
fed_id NUMBER(3) not null,
tt_id  NUMBER(3) not null,
gol_id  NUMBER(2) not null,
gol_current NUMBER(10),
pll_current NUMBER(10),
status    NUMBER(10)
);

ALTER TABLE GOL_BIAS_CURRENT_DAT ADD CONSTRAINT GOL_BIAS_CURRENT_DAT_pk PRIMARY KEY (rec_id, fed_id, tt_id, gol_id );



/* FE DAQ */


CREATE TABLE FE_DAQ_CONFIG (
config_id 	NUMBER(10) NOT NULL,
tag VARCHAR2(20) not null,
version NUMBER(10) not null,
ped_id NUMBER(10), 
del_id NUMBER(10), 
wei_id NUMBER(10),
bxt_id NUMBER(10),
btt_id NUMBER(10),
 TR_BXT_ID                                          NUMBER(10),
 TR_BTT_ID                                          NUMBER(10),
 TBY_ID                                             NUMBER(10),
 VFE_ID                                             NUMBER(10),
 GOL_ID                                             NUMBER(10),
USER_COMMENT VARCHAR2(100),
db_timestamp          TIMESTAMP DEFAULT SYSTIMESTAMP NOT NULL
);

ALTER TABLE FE_DAQ_CONFIG ADD CONSTRAINT FE_DAQ_CONFIG_pk PRIMARY KEY (config_id);
ALTER TABLE FE_DAQ_CONFIG ADD CONSTRAINT FE_DAQ_CONFIG_uk UNIQUE(tag,version);
ALTER TABLE FE_DAQ_CONFIG ADD CONSTRAINT FE_DAQ_CONFIG_fk1 FOREIGN KEY (ped_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE FE_DAQ_CONFIG ADD CONSTRAINT FE_DAQ_CONFIG_fk2 FOREIGN KEY (del_ID) REFERENCES COND2CONF_INFO (REC_ID);
ALTER TABLE FE_DAQ_CONFIG ADD CONSTRAINT FE_DAQ_CONFIG_fk3 FOREIGN KEY (wei_ID) REFERENCES COND2CONF_INFO (REC_ID);

CREATE SEQUENCE FE_DAQ_CONDFIG_SQ INCREMENT BY 1 START WITH 1;



---- Some mapping tables

create table ECAL_FED_DEF(
	DEF_ID NUMBER NOT NULL,
	HOST VARCHAR2(100) NOT NULL,
	SLOT NUMBER NOT NULL,
	BOARD_ID NUMBER, 
	FED_ID NUMBER NOT NULL
);

-- DOESN'T WORK FOR PRESHOWER
ALTER TABLE ECAL_FED_DEF ADD CONSTRAINT ecal_fed_def_pk  PRIMARY KEY (DEF_ID);
-- ALTER TABLE ECAL_FED_DEF ADD CONSTRAINT ecal_fed_def_uk  UNIQUE (host,slot);

-- NOT USED ANYMORE
CREATE SEQUENCE ecal_fed_def_sq INCREMENT BY 1 START WITH 1;
CREATE trigger ecal_fed_def_trg
before insert on ECAL_FED_DEF
for each row
begin
select ecal_fed_def_sq.NextVal into :new.def_id from dual;
end;
/

create table ECAL_FED_TO_SUPERMODULE(
        DEF_ID          NUMBER NOT NULL,
        FED_ID          NUMBER NOT NULL,
        CONSTR_ID       NUMBER NOT NULL,
        PSEUDO_SLOT_ID  NUMBER NOT NULL,
        GEOM_ID         VARCHAR(10)
);

ALTER TABLE ECAL_FED_TO_SUPERMODULE ADD CONSTRAINT ECAL_FED_TO_SUPERMODULE_PK  PRIMARY KEY (DEF_ID);
ALTER TABLE ECAL_FED_TO_SUPERMODULE ADD CONSTRAINT ECAL_FED_TO_SUPERMODULE_UK1 UNIQUE (FED_ID);
ALTER TABLE ECAL_FED_TO_SUPERMODULE ADD CONSTRAINT ECAL_FED_TO_SUPERMODULE_UK2 UNIQUE (CONSTR_ID);
ALTER TABLE ECAL_FED_TO_SUPERMODULE ADD CONSTRAINT ECAL_FED_TO_SUPERMODULE_UK3 UNIQUE (PSEUDO_SLOT_ID);

CREATE SEQUENCE ECAL_FED_TO_SUPERMODULE_SQ INCREMENT BY 1 START WITH 1;
CREATE trigger ECAL_FED_TO_SUPERMODULE_TRG
before insert on ECAL_FED_TO_SUPERMODULE
for each row
begin
select ECAL_FED_TO_SUPERMODULE_SQ.NextVal into :new.def_id from dual;
end;
/
CREATE OR REPLACE procedure cond2conf_autoinsert
(rec_type in varchar2, rec_id in NUMBER, tag in varchar2, version in number)
IS

 sql_str VARCHAR(1000);
 location_id number(10); 
 type_id number;
  loca varchar2(10);
  short_descr varchar2(30) ;

BEGIN

  sql_str := 'SELECT def_id from cond2conf_type_def where rec_type=:1 ';
    EXECUTE IMMEDIATE sql_str INTO type_id using rec_type ;
  loca:='P5';
  sql_str := 'SELECT def_id from location_def where location=:1 ';
    EXECUTE IMMEDIATE sql_str INTO location_id using loca ;
  short_descr:= tag || '_' || cast (version as varchar2) ;
  sql_str :='Insert into COND2CONF_INFO (rec_id,REC_TYPE_ID, LOCATION_ID, short_desc ) values (:1, :2, :3, :4) ';
    EXECUTE IMMEDIATE sql_str  using rec_id, type_id, location_id, short_descr  ;
end;
/
show errors;




create OR REPLACE TRIGGER cond2conf_auto_tg
  BEFORE INSERT ON delays_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('DELAY_OFFSETS', :newiov.rec_id, :newiov.tag, :newiov.version)
/
create OR REPLACE TRIGGER cond2conf_auto2_tg
  BEFORE INSERT ON pedestal_offsets_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('PEDESTAL_OFFSETS', :newiov.rec_id, :newiov.tag, :newiov.version)
/
create OR REPLACE TRIGGER cond2conf_auto_tg3
  BEFORE INSERT ON weights_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('DCC_WEIGHTS', :newiov.rec_id, :newiov.tag, :newiov.version)
/
create OR REPLACE TRIGGER cond2conf_auto_tg4
  BEFORE INSERT ON bad_crystals_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('BAD_CRYSTALS', :newiov.rec_id, :newiov.tag, :newiov.version)
/
create OR REPLACE TRIGGER cond2conf_auto_tg5
  BEFORE INSERT ON bad_tt_info
  REFERENCING NEW AS newiov
  FOR EACH ROW
  CALL cond2conf_autoinsert('BAD_TT', :newiov.rec_id, :newiov.tag, :newiov.version)
/