-- Not tested script for population of testing entries to ovini database

CREATE DATABASE ovini OWNER www-data

--
-- PostgreSQL database dump
--

SET client_encoding = 'UTF8';
SET standard_conforming_strings = off;
SET check_function_bodies = false;
SET client_min_messages = warning;
SET escape_string_warning = off;

SET search_path = public, pg_catalog;

ALTER TABLE ONLY public.jobs DROP CONSTRAINT jobs_pkey;
DROP TABLE public.jobs;
DROP SCHEMA public;
--
-- Name: public; Type: SCHEMA; Schema: -; Owner: postgres
--

CREATE SCHEMA public;


ALTER SCHEMA public OWNER TO postgres;

--
-- Name: SCHEMA public; Type: COMMENT; Schema: -; Owner: postgres
--

COMMENT ON SCHEMA public IS 'standard public schema';


SET search_path = public, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: jobs; Type: TABLE; Schema: public; Owner: postgres; Tablespace:
--

CREATE TABLE jobs (
    id character varying(10) DEFAULT 0 NOT NULL,
    type character varying(100),
    status character varying(100),
    created character varying(100),
    config text
);


ALTER TABLE public.jobs OWNER TO postgres;

--
-- Data for Name: jobs; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY jobs (id, type, status, created, config) FROM stdin;
2	phonon	finished	2009-06-11 16:47:52.790506	phonons of Ni at gamma\n &inputph\n  tr2_ph=1.0d-16,\n  prefix='ni',\n  ldisp=.true.,\n  nq1=2,\n  nq2=2,\n  nq3=2,\n  amass(1)=58.6934,\n  outdir='temp/',\n  fildyn='ni.dyn',\n /\n
1	electron	finished	2009-06-11 13:35:34.081567	 &control\n    calculation='scf'\n    restart_mode='from_scratch',\n    tprnfor = .true.\n    prefix='ni',\n    pseudo_dir = './',\n    outdir='temp/'\n /\n &system    \n    ibrav=2, \n    celldm(1) =6.65, \n    nat=  1, \n    ntyp= 1,\n    nspin=2,\n    starting_magnetization(1)=0.5,\n    degauss=0.02,\n    smearing='gauss',\n    occupations='smearing',\n    ecutwfc =27.0\n    ecutrho =300.0\n /\n &electrons\n    conv_thr =  1.0d-8\n    mixing_beta = 0.7\n /\nATOMIC_SPECIES\n Ni  26.98  Ni.pbe-nd-rrkjus.UPF\nATOMIC_POSITIONS\n Ni 0.00 0.00 0.00 \nK_POINTS AUTOMATIC\n4 4 4 1 1 1\n
3	phonon	running	2009-07-16 23:06:15.024084	phonons of Ni at gamma\n &inputph\n  tr2_ph=1.0d-16,\n  prefix='ni',\n  ldisp=.true.,\n  nq1=2,\n  nq2=2,\n  nq3=2,\n  amass(1)=58.6934,\n  outdir='temp/',\n  fildyn='ni.dyn',\n /\n
\.


--
-- Name: jobs_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace:
--

ALTER TABLE ONLY jobs
    ADD CONSTRAINT jobs_pkey PRIMARY KEY (id);


--
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

