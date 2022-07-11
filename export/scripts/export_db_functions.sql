
-- copies newly created NOAA active regions from private production series to
-- public read-only series used by the export system (processing.elements.js)
CREATE OR REPLACE FUNCTION jsoc.insert_noaa_active_region() RETURNS TRIGGER AS
$$
BEGIN
  IF (TG_OP='INSERT' AND NEW.observationtime > 0 AND NEW.regionnumber > 0) THEN
      INSERT INTO jsoc.noaa_active_regions(recnum, sunum, slotnum, sessionid, sessionns, observationtime, regionnumber, zurichclass, magnetictype, spotcount, area, latitudehg, longitudehg, longitudecm, longitudinalextent) VALUES (NEW.recnum, NEW.sunum, NEW.slotnum, 0, 'NA', NEW.observationtime, NEW.regionnumber, NEW.zurichclass, NEW.magnetictype, NEW.spotcount, NEW.area, NEW.latitudehg, NEW.longitudehg, NEW.longitudecm, NEW.longitudinalextent);
  END IF;
RETURN NEW;
END
$$
LANGUAGE plpgsql;

-- db user rick needs to be able to run jsoc.insert_noaa_active_region
GRANT USAGE ON SCHEMA jsoc TO rick;
-- db user rick needs to be able to insert into jsoc.noaa_active_regions
GRANT INSERT ON jsoc.noaa_active_regions TO rick;
-- create trigger function called when a row is inserted into su_rsb.noaa_activeregions
-- trigger function copies the inserted row into jsoc.noaa_active_regions
DROP TRIGGER IF EXISTS insert_noaa_active_region ON su_rsb.noaa_activeregions;
CREATE TRIGGER insert_noaa_active_region AFTER INSERT ON su_rsb.noaa_activeregions FOR EACH ROW EXECUTE PROCEDURE jsoc.insert_noaa_active_region()
