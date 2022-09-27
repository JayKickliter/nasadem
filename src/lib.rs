//! Parsers for NASA Digital Elevation Model.

use byteorder::{BigEndian as BE, ReadBytesExt};
use geo_types::{LineString, Point, Polygon};
use std::io::{Error as IoError, Read};

type DEMMatrix<T> = Vec<T>;

#[derive(Debug)]
pub struct NASADEM {
    southwest_corner: Point<i32>,
    elevation: Option<DEMMatrix<u16>>,
    water: Option<DEMMatrix<bool>>,
}

impl NASADEM {
    pub fn new(southwest_corner: Point<i32>) -> Self {
        Self {
            southwest_corner,
            elevation: None,
            water: None,
        }
    }

    pub fn add_elevation(&mut self, mut src: impl Read) -> Result<&mut Self, IoError> {
        let mut elev_samples = Vec::with_capacity(3601 * 3601);
        let mut idx = 0_usize;
        for y in (0..3601).rev() {
            let lat_b = self.southwest_corner.y() as f64 + y as f64 / 3601.0;
            let lat_t = lat_b + 1.0 / 3601.0;
            debug_assert!(lat_t <= (self.southwest_corner.y() + 1) as f64);
            debug_assert!(lat_b < (self.southwest_corner.y() + 1) as f64);
            debug_assert!(lat_b >= self.southwest_corner.y() as f64);
            debug_assert!(lat_t > self.southwest_corner.y() as f64);
            for x in 0..3601 {
                let lon_l = self.southwest_corner.x() as f64 + x as f64 / 3601.0;
                let lon_r = lon_l + 1.0 / 3601.0;
                debug_assert!(lon_r <= (self.southwest_corner.x() + 1) as f64);
                debug_assert!(lon_l < (self.southwest_corner.x() + 1) as f64);
                debug_assert!(lon_r > self.southwest_corner.x() as f64);
                debug_assert!(lon_l >= self.southwest_corner.x() as f64);
                let sample = src.read_u16::<BE>()?;
                elev_samples.push(sample);
                debug_assert_eq!(
                    (idx, idx_to_pont(&self.southwest_corner, idx)),
                    (idx, Point::new(lon_l, lat_b))
                );
                idx += 1;
            }
        }
        debug_assert_eq!(elev_samples.len(), 3601 * 3601);
        self.elevation = Some(elev_samples);
        Ok(self)
    }

    pub fn add_water(&mut self, mut src: impl Read) -> Result<&mut Self, IoError> {
        let mut water_samples = Vec::with_capacity(3601 * 3601);
        for _i in 0..3601 {
            for _j in 0..3601 {
                let sample = src.read_u8()?;
                debug_assert!(sample == 0 || sample == 255);
                water_samples.push(sample == 255);
            }
        }
        debug_assert_eq!(water_samples.len(), 3601 * 3601);
        self.water = Some(water_samples);
        Ok(self)
    }

    pub fn iter(&'_ self) -> impl Iterator<Item = DEMBox> + '_ {
        Iter { dem: self, idx: 0 }
    }
}

pub fn idx_to_pont(sw_corner: &Point<i32>, idx: usize) -> Point<f64> {
    debug_assert!(idx < 3601 * 3601);
    let y = 3600 - (idx / 3601);
    let lat_south = sw_corner.y() as f64 + y as f64 / 3601.0;
    let x = idx % 3601;
    let lon_west = sw_corner.x() as f64 + x as f64 / 3601.0;
    Point::new(lon_west, lat_south)
}

struct Iter<'a> {
    dem: &'a NASADEM,
    idx: usize,
}

impl<'a> Iterator for Iter<'a> {
    type Item = DEMBox;

    fn next(&mut self) -> Option<DEMBox> {
        if self.idx < 3601 * 3601 {
            let southwest_corner = idx_to_pont(&self.dem.southwest_corner, self.idx);
            let elevation = self.dem.elevation.as_ref().map(|e| e[self.idx]);
            let is_water = self.dem.water.as_ref().map(|w| w[self.idx]);
            self.idx += 1;
            Some(DEMBox {
                southwest_corner,
                elevation,
                is_water,
            })
        } else {
            None
        }
    }
}

pub struct DEMBox {
    southwest_corner: Point<f64>,
    elevation: Option<u16>,
    is_water: Option<bool>,
}

impl DEMBox {
    pub fn polygon(&self) -> Polygon {
        let lat_south = self.southwest_corner.y();
        let lat_north = lat_south + 1.0 / 3601.0;
        let lon_west = self.southwest_corner.x();
        let lon_east = lon_west + (1.0 / 3601.0);
        Polygon::new(
            LineString::from(vec![
                (lon_west, lat_south),
                (lon_east, lat_south),
                (lon_east, lat_north),
                (lon_west, lat_north),
                (lon_west, lat_south),
            ]),
            Vec::new(),
        )
    }

    pub fn southwest_corner(&self) -> &Point {
        &self.southwest_corner
    }

    pub fn elevation(&self) -> Option<u16> {
        self.elevation
    }

    pub fn is_water(&self) -> Option<bool> {
        self.is_water
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hextree::{compaction::EqCompactor, h3ron, HexTreeMap};
    use std::{
        fs::File,
        io::{BufReader, BufWriter},
    };

    #[test]
    fn test_new() {
        let elevation_src = BufReader::new(
            File::open(format!(
                "{}/local/NASADEM_HGT_n38w106/n38w106.hgt",
                std::env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        );

        let water_src = BufReader::new(
            File::open(format!(
                "{}/local/NASADEM_HGT_n38w106/n38w106.swb",
                std::env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        );

        let mut dem = NASADEM::new(Point::new(-106, 38));
        dem.add_elevation(elevation_src).unwrap();
        dem.add_water(water_src).unwrap();

        let mut iter = dem.iter();
        let dbox_0_0 = iter.next().unwrap();
        assert_eq!(
            dbox_0_0.southwest_corner(),
            &Point::new(-106.0, 38.99972229936129)
        );
    }

    #[test]
    fn test_hex_map() {
        let elevation_src = BufReader::new(
            File::open(format!(
                "{}/local/NASADEM_HGT_n38w106/n38w106.hgt",
                std::env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        );

        let mut dem = NASADEM::new(Point::new(-106, 38));
        dem.add_elevation(elevation_src).unwrap();

        let mut elev_map = HexTreeMap::with_compactor(EqCompactor);

        let mut pre_compaction_cell_count = 0;
        for (n, dem_box) in dem.iter().enumerate() {
            let elev = dem_box.elevation().unwrap();
            for cell in &h3ron::polygon_to_cells(&dem_box.polygon(), 14).unwrap() {
                elev_map.insert(cell, elev);
                pre_compaction_cell_count += 1;
            }
            if n > 0 && n % 3601 == 0 {
                println!(
                    "sample {} ({:.02}%), total cells {}",
                    n,
                    (n * 100) as f64 / (3601.0 * 3601.0),
                    pre_compaction_cell_count
                );
            }
        }

        let out = BufWriter::new(
            File::create(format!(
                "{}/local/NASADEM_HGT_n38w106/n38w106.res14.bincode.hexmap",
                std::env!("CARGO_MANIFEST_DIR")
            ))
            .unwrap(),
        );

        bincode::serialize_into(out, &elev_map).unwrap();

        println!(
            "map len {}, total cells {}",
            elev_map.len(),
            pre_compaction_cell_count
        );
        assert!(elev_map.len() < pre_compaction_cell_count);
    }
}
