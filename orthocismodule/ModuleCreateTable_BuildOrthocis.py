import sys
try:
        import MySQLdb
        bHaveMySQLdb=True
except ImportError:
        bHaveMySQLdb=False

from ModuleFonction_BuildOrthocis import *

def createComConnexTable(cctablenm,ccgid,optsql,host,database,user,pw,db,sVers,bHavedb):

                sSQLcreateComponentConnexTable="drop table if exists "+cctablenm+";"
                sSQLcreateComponentConnexTable+="create table "+cctablenm+" "
                sSQLcreateComponentConnexTable+="("+ccgid+" bigint(14), "
                sSQLcreateComponentConnexTable+="component_id bigint(14) default NULL,primary key("+ccgid+"));"
                sSQLcreateComponentConnexTable+="alter table "+cctablenm+ " ADD INDEX idx_A  (component_id );"
                #print sSQLcreateComponentConnexTable
                SubmitMySQLCommand(sSQLcreateComponentConnexTable,host,database,user,pw,db,bHaveMySQLdb,False)

                sSQLpopulateComponentConnexTable="insert IGNORE into "+cctablenm+" ("+ccgid+") select "+ccgid+" from "+sVers+"_gene where "+optsql+" ;"
                #print sSQLpopulateComponentConnexTable 
                SubmitMySQLCommand(sSQLpopulateComponentConnexTable,host,database,user,pw,db,bHavedb,False)

		


#Content of the tool table:
def GetContentOfToolTable():
	return {
	"TESS":(1,"score"),
	"LOGOL":(2,"undef"),
	"RSAT":(3,"Pvalue")}

#Prepare SQL syntax to create table specific to species
def DefineSQLRequestSpecificToSpecie(stringVersion,strgid,upstreamsize):
	#print stringVersion,strgid
	return [
"drop table if exists "+stringVersion+"_scan_"+strgid+""";
CREATE TABLE """+stringVersion+"_scan_"+strgid+""" (
  scan_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  motif_id int(10) unsigned NOT NULL,
  gene_id bigint(14) unsigned NOT NULL,
  \`match\` tinyint(1) DEFAULT 1,
  job_id int(10) unsigned NOT NULL,
  start_match int(10),
  end_match int(10),
  sequence_match varchar(20),
  strand_match tinyint(1),
  score float(10) DEFAULT NULL,
  Pvalue double DEFAULT NULL,
  region int(3),
  PRIMARY KEY (scan_id),
  KEY (motif_id),
  KEY (gene_id),
  KEY (job_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
 "drop table if exists "+stringVersion+"_sequence_"+strgid+""";
  CREATE TABLE """+stringVersion+"_sequence_"+strgid+""" (
  gene_id bigint(14) unsigned NOT NULL,
  sequence MEDIUMTEXT,
  upstreamsize bigint(9), 
  genesize bigint(9),
  genefragsize bigint(9),
  isnegativefragsize int(10) DEFAULT 0,	
  PRIMARY KEY (gene_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
 ,
  "drop table if exists "+stringVersion+"_transcript_"+strgid+""";
  CREATE TABLE """+stringVersion+"_transcript_"+strgid+"""(
  transcript_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  gene_id bigint(14) unsigned DEFAULT NULL,
  analysis_id smallint(5) unsigned NOT NULL,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(10) unsigned NOT NULL,
  seq_region_end int(10) unsigned NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  display_xref_id int(10) unsigned DEFAULT NULL,
  src_fixme varchar(40) NOT NULL,
  biotype varchar(40) NOT NULL,
  description text,
  is_current tinyint(1) NOT NULL DEFAULT '1',
  canonical_translation_id int(10) unsigned DEFAULT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) unsigned NOT NULL DEFAULT '1',
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  PRIMARY KEY (transcript_id),
  UNIQUE KEY canonical_translation_idx_"""+strgid+""" (canonical_translation_id),
  KEY seq_region_idx (seq_region_id,seq_region_start),
  KEY gene_index (gene_id),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id),
  KEY stable_id_idx (stable_id,version)
) ENGINE=MyISAM AUTO_INCREMENT=26741 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"_assembly_"+strgid+""";
  CREATE TABLE """+stringVersion+"_assembly_"+strgid+""" (
  asm_seq_region_id int(10) unsigned NOT NULL,
  cmp_seq_region_id int(10) unsigned NOT NULL,
  asm_start int(10) NOT NULL,
  asm_end int(10) NOT NULL,
  cmp_start int(10) NOT NULL,
  cmp_end int(10) NOT NULL,
  ori tinyint(4) NOT NULL,
  UNIQUE KEY all_idx_"""+strgid+""" (asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori),
  KEY cmp_seq_region_idx (cmp_seq_region_id),
  KEY asm_seq_region_idx (asm_seq_region_id,asm_start)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"_coord_system_"+strgid+""";
  CREATE TABLE """+stringVersion+"_coord_system_"+strgid+""" (
  coord_system_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  species_id int(10) unsigned NOT NULL DEFAULT '1',
  name varchar(40) NOT NULL,
  version varchar(255) DEFAULT NULL,
  rank int(11) NOT NULL,
  attrib set('default_version','sequence_level') DEFAULT NULL,
  PRIMARY KEY (coord_system_id),
  UNIQUE KEY rank_idx_"""+strgid+""" (rank,species_id),
  UNIQUE KEY name_idx_"""+strgid+""" (name,version,species_id),
  KEY species_idx (species_id)
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"_dna_"+strgid+""";
  CREATE TABLE """+stringVersion+"_dna_"+strgid+""" (
  seq_region_id int(10) unsigned NOT NULL,
  sequence longtext NOT NULL,
  PRIMARY KEY (seq_region_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=750000 AVG_ROW_LENGTH=19000;"""
,
"drop table if exists "+stringVersion+"_gene_"+strgid+""";
  CREATE TABLE """+stringVersion+"_gene_"+strgid+""" (
  gene_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  biotype varchar(40) NOT NULL,
  analysis_id smallint(5) unsigned NOT NULL,
  seq_region_id int(10) unsigned NOT NULL,
  seq_region_start int(10) unsigned NOT NULL,
  seq_region_end int(10) unsigned NOT NULL,
  seq_region_strand tinyint(2) NOT NULL,
  display_xref_id int(10) unsigned DEFAULT NULL,
  source varchar(20) NOT NULL,
  description text,
  is_current tinyint(1) NOT NULL DEFAULT '1',
  canonical_transcript_id int(10) unsigned NOT NULL,
  stable_id varchar(128) DEFAULT NULL,
  version smallint(5) unsigned NOT NULL DEFAULT '1',
  created_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  modified_date datetime NOT NULL DEFAULT '0000-00-00 00:00:00',
  PRIMARY KEY (gene_id),
  KEY seq_region_idx (seq_region_id,seq_region_start),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id),
  KEY stable_id_idx (stable_id,version),
  KEY canonical_transcript_id_idx (canonical_transcript_id)
) ENGINE=MyISAM AUTO_INCREMENT=21735 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"_seq_region_"+strgid+""";
  CREATE TABLE """+stringVersion+"_seq_region_"+strgid+""" (
  seq_region_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  name varchar(40) NOT NULL,
  coord_system_id int(10) unsigned NOT NULL,
  length int(10) unsigned NOT NULL,
  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_cs_idx_"""+strgid+""" (name,coord_system_id),
  KEY name_idx (name),
  KEY cs_idx (coord_system_id)
) ENGINE=MyISAM AUTO_INCREMENT=102843 DEFAULT CHARSET=latin1;"""
]

#Prepare SQL syntax to create generic table
def DefineSQLRequestGeneric(stringVersion,upstreamsize):
	return [
"drop table if exists "+stringVersion+"""_genome_db;
  CREATE TABLE """+stringVersion+"""_genome_db (
  genome_db_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  taxon_id int(10) unsigned DEFAULT NULL,
  name varchar(40) NOT NULL DEFAULT '',
  assembly varchar(100) NOT NULL DEFAULT '',
  assembly_default tinyint(1) DEFAULT '1',
  genebuild varchar(100) NOT NULL DEFAULT '',
  locator varchar(400) DEFAULT NULL,
  PRIMARY KEY (genome_db_id),
  UNIQUE KEY name (name,assembly,genebuild)
) ENGINE=MyISAM AUTO_INCREMENT=144 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_HGNC;
  CREATE TABLE """+stringVersion+"""_HGNC (
  HGNC_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  fullName varchar(100) NOT NULL,
  shortName varchar(100) NOT NULL,
  PRIMARY KEY (HGNC_id)
  )ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_asso_HGNC_clique;
CREATE TABLE """+stringVersion+"""_asso_HGNC_clique (
  HGNC_id int(10) unsigned NOT NULL,
  clique_id int(10) unsigned NOT NULL,
  UNIQUE KEY asso_hgnc_gene1 (HGNC_id,clique_id)
  )ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_job;
CREATE TABLE """+stringVersion+"""_job (
  job_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  tool_id int(10) unsigned DEFAULT NULL,
  motif_id int(10) unsigned DEFAULT NULL,
  param varchar(250) DEFAULT NULL,
  PRIMARY KEY (job_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_parameter;
CREATE TABLE """+stringVersion+"""_parameter (
  parameter_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  tool_id int(10) unsigned DEFAULT NULL,
  name varchar(40) NOT NULL DEFAULT '',
  \`desc\` varchar(100),
  PRIMARY KEY (parameter_id)
)ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_tool;
CREATE TABLE """+stringVersion+"""_tool (
  tool_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  name varchar(36) NOT NULL,
  Tag varchar(50),
  PRIMARY KEY (tool_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_motif;
CREATE TABLE """+stringVersion+"""_motif (
  motif_id int(10) NOT NULL AUTO_INCREMENT,
  name varchar(40) NOT NULL DEFAULT '',
  content varchar(10000) DEFAULT NULL,
  score float DEFAULT NULL,
  size int(10),
  PRIMARY KEY (motif_id),
  UNIQUE KEY kn (name)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_homolink;
CREATE TABLE """+stringVersion+"""_homolink ( 
  first bigint(14) unsigned NOT NULL, 
  second bigint(14) unsigned NOT NULL, 
  UNIQUE KEY link (first,second) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_globalhomolinkstart;
CREATE TABLE """+stringVersion+"""_globalhomolinkstart ( 
  first bigint(14) unsigned NOT NULL, 
  second bigint(14) unsigned NOT NULL, 
  UNIQUE KEY link (first,second) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,

"drop table if exists "+stringVersion+"""_global_homolink;
CREATE TABLE """+stringVersion+"""_global_homolink ( 
  first bigint(14) unsigned NOT NULL, 
  second bigint(14) unsigned NOT NULL, 
  UNIQUE KEY link (first,second) 
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_gene;
CREATE TABLE """+stringVersion+"""_gene (
  gene_id bigint(14) unsigned NOT NULL AUTO_INCREMENT,
  stable_id varchar(128) NOT NULL,
  genome_db_id int(14) unsigned DEFAULT NULL,
  upstream_size int("""+str(len(str(upstreamsize)))+""") DEFAULT NULL,
  global_id bigint(14) unsigned NOT NULL,
  ismany int(1) DEFAULT 0,
  PRIMARY KEY (gene_id),
  KEY (stable_id),
  KEY (global_id) 
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_homology;
CREATE TABLE """+stringVersion+"""_homology (
  homology_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  method_link_species_set_id int(10) unsigned NOT NULL,
  description enum('ortholog_one2one','apparent_ortholog_one2one','ortholog_one2many','ortholog_many2many','within_species_paralog','other_paralog','putative_gene_split','contiguous_gene_split','between_species_paralog','possible_ortholog','UBRH','BRH','MBRH','RHS','projection_unchanged','projection_altered')DEFAULT NULL,
  PRIMARY KEY (homology_id),
  KEY method_link_species_set_idx (method_link_species_set_id),
  KEY descriptionidx (description)
) ENGINE=MyISAM AUTO_INCREMENT=400000066 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_homology_full;
CREATE TABLE """+stringVersion+"""_homology_full (
  homology_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  method_link_species_set_id int(10) unsigned NOT NULL,
  description enum('ortholog_one2one','apparent_ortholog_one2one','ortholog_one2many','ortholog_many2many','within_species_paralog','other_paralog','putative_gene_split','contiguous_gene_split','between_species_paralog','possible_ortholog','UBRH','BRH','MBRH','RHS','projection_unchanged','projection_altered')DEFAULT NULL,
  subtype varchar(40) NOT NULL DEFAULT '',
  dn float(10,5) DEFAULT NULL,
  ds float(10,5) DEFAULT NULL,
  n float(10,1) DEFAULT NULL,
  s float(10,1) DEFAULT NULL,
  lnl float(10,3) DEFAULT NULL,
  threshold_on_ds float(10,5) DEFAULT NULL,
  ancestor_node_id int(10) unsigned NOT NULL,
  tree_node_id int(10) unsigned NOT NULL,
  PRIMARY KEY (homology_id)
) ENGINE=MyISAM AUTO_INCREMENT=400000066 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_homology_member;
CREATE TABLE """+stringVersion+"""_homology_member (
  homology_id int(10) unsigned NOT NULL,
  member_id int(10) unsigned NOT NULL,
  PRIMARY KEY (homology_id,member_id),
  KEY homology_id (homology_id),
  KEY member_id (member_id)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_member;
CREATE TABLE """+stringVersion+"""_member (
  member_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  stable_id varchar(128) NOT NULL,
  genome_db_id int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (member_id),
  KEY stable_id (stable_id),
  KEY genome_db_id (genome_db_id)
) ENGINE=MyISAM AUTO_INCREMENT=400000065 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_method_link_species_set;
CREATE TABLE """+stringVersion+"""_method_link_species_set (
  method_link_species_set_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  method_link_id int(10) unsigned DEFAULT NULL,
  species_set_id int(10) unsigned NOT NULL DEFAULT '0',
  name varchar(255) NOT NULL DEFAULT '',
  source varchar(255) NOT NULL DEFAULT 'ensembl',
  url varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (method_link_species_set_id),
  UNIQUE KEY method_link_id (method_link_id,species_set_id)
) ENGINE=MyISAM AUTO_INCREMENT=50047 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_species_set;
CREATE TABLE """+stringVersion+"""_species_set (
  species_set_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  genome_db_id int(10) unsigned DEFAULT NULL,
  UNIQUE KEY species_set_id (species_set_id,genome_db_id),
  KEY genome_db_id (genome_db_id)
) ENGINE=MyISAM AUTO_INCREMENT=35266 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_method_link;
CREATE TABLE """+stringVersion+"""_method_link (
  method_link_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  type varchar(50) NOT NULL DEFAULT '',
  class varchar(50) NOT NULL DEFAULT '',
  PRIMARY KEY (method_link_id),
  UNIQUE KEY type (type)
) ENGINE=MyISAM AUTO_INCREMENT=503 DEFAULT CHARSET=latin1;"""
,
"drop table if exists "+stringVersion+"""_method_link;
CREATE TABLE """+stringVersion+"""_method_link (
  method_link_id int(10) unsigned NOT NULL AUTO_INCREMENT,
  type varchar(50) NOT NULL DEFAULT '',
  class varchar(50) NOT NULL DEFAULT '',
  PRIMARY KEY (method_link_id),
  UNIQUE KEY type (type)
) ENGINE=MyISAM AUTO_INCREMENT=503 DEFAULT CHARSET=latin1;"""
]
