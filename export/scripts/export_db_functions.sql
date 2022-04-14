
-- copies newly created NOAA active regions from private production series to
-- public read-only series used by the export system (processing.elements.js)
CREATE OR REPLACE FUNCTION jsoc.insert_noaa_active_region() RETURNS TRIGGER AS
$$
BEGIN
  IF (TG_OP='INSERT' AND NEW.observationtime > 0 AND NEW.regionnumber > 0) THEN
      INSERT INTO jsoc.noaa_active_regions(observationtime, regionnumber, zurichclass, magnetictype, spotcount, area, latitudehg, longitudehg, longitudecm, longitudinalextent) VALUES (NEW.observationtime, NEW.regionnumber, NEW.zurichclass, NEW.magnetictype, NEW.spotcount, NEW.area, NEW.latitudehg, NEW.longitudehg, NEW.longitudecm, NEW.longitudinalextent);
  END IF;
RETURN NEW;
END
$$
LANGUAGE plpgsql;
