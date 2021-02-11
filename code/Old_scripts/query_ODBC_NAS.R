ch<-odbcConnect("PostgreSQL30")  # définit le pilote ODBC

my.sp.rel<-sqlQuery(ch,paste("
                             SELECT public.observations_stable.obs_id, public.observations_stable.locality_id, public.observations_stable.project_id, 
                             public.observations_stable.input_id, public.observations_stable.releve_type, public.observations_stable.obs_date, 
                             public.observations_stable.date_precision, public.observations_stable.date_expert, public.observations_stable.introduced, 
                             public.observations_stable.introduced_expert, public.observations_stable.x_swiss, public.observations_stable.y_swiss, 
                             public.observations_stable.xy_type, public.observations_stable.xy_precision, public.observations_stable.count_unit,
                             public.observations_stable.abundance_code, public.observations_stable.v_taxon, public.observations_stable.v_accepted_taxon_id, 
                             public.observations_stable.v_validation_status, public.observations_stable.v_observers, public.observations_stable.v_native_status, 
                             public.observations_stable.v_presence_status, public.observations_stable.v_introduction_status, public.observations_stable.v_doubt_status, 
                             public.observations_stable.v_doubt_status, public.observations_stable.v_xy_radius, ST_AsText(ST_SimplifyPreserveTopology(ST_TRANSFORM(public.observations_stable.geom_storage, 21781),50)),
                             ST_AsText(ST_SimplifyPreserveTopology(ST_TRANSFORM(public.observations_stable.v_geom_buffer, 21781),50))
                             FROM public.observations_stable 
                             WHERE (((public.observations_stable.v_accepted_taxon_id)=",sp.rel[i,2],")  AND (public.observations_stable.presence = 681));
                             ",sep=""))  #extract the database for the given species id_number
my.sp<-rbind(my.sp,my.sp.rel)

odbcCloseAll()

####################################################

x<-sp.no #code of the taxon name
ch<-odbcConnect("PostgreSQL35W") # Extraction of the database for a given number (corresponding to the code of the taxon name)

sp.ws<-sqlQuery(ch,paste("
                         SELECT public.welten_sutter.no_surface, public.welten_sutter.no_welten_sutter, public.welten_sutter.no_nom,
                         public.welten_sutter.code    
                         FROM public.welten_sutter 
                         WHERE (((public.welten_sutter.no_nom)=",x,")  AND (public.welten_sutter.code < 6));
                         ",sep=""))  #extract the WS sectors for the given taxon id_number

sp.indigenat<-sqlQuery(ch,paste("
                                SELECT public.nom_indigenat_biogeo.no_nom, public.nom_indigenat_biogeo.no_region, 
                                public.nom_indigenat_biogeo.indigenat_biogeo
                                FROM public.nom_indigenat_biogeo 
                                WHERE (((public.nom_indigenat_biogeo.no_nom)=",x,")  AND (public.nom_indigenat_biogeo.indigenat_biogeo < 513));
                                ",sep=""))  #extract bioregions in which taxon is indigeneous


odbcCloseAll()
