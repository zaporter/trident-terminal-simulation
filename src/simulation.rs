use crate::*;

use kiss3d::camera::FirstPerson;
use kiss3d::nalgebra::{Vector3, UnitQuaternion, Point3};
use kiss3d::ncollide3d::math::Translation;
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use kiss3d::light::Light;

use stopwatch::Stopwatch;

#[derive(Debug, Clone)]
pub struct SimulationConfig{
    pub g : SimulationFloat,
    pub do_visuals : bool,
}

pub trait Simulatable {
    fn tick(&mut self, dt_s : &SimulationFloat, gravity_applying_masses : &Vec<&PointMass>);
    fn calculate_momentum(&self) -> Vec3;
    fn calculate_total_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) -> SimulationFloat;
    fn update_scene_node_position(&mut self);
    fn draw_historical_positions(&self, window : &mut Window);
}

#[derive(Default, Debug, Clone)]
pub struct Simulation {
    pub point_masses : Vec<PointMass>,
    pub ropes : Vec<Rope>,
    pub simulation_elapsed_time_s : SimulationFloat,
}


impl Simulation {

    pub fn run_till<F>(&mut self, simulation_jiffy_s : &SimulationFloat, criteria : F) where
        F: Fn(SimulationFloat)->bool {
        while criteria(self.simulation_elapsed_time_s.clone()) {

            let masses_copy = self.point_masses.clone();
            let mut grav_masses = Vec::new();
            for mass in &masses_copy {
                if mass.applies_gravity{
                    grav_masses.push(mass);
                }
            }
            for point in &mut self.point_masses {
                point.tick(&simulation_jiffy_s, &grav_masses);
            }
            self.simulation_elapsed_time_s += simulation_jiffy_s.clone();
        }
    }
    pub fn calculate_total_energy_j(&self) -> SimulationFloat{
        let mut total_energy_j = SimulationFloat::default();
        let masses_copy = self.point_masses.clone();
        let mut grav_masses = Vec::new();
        for mass in &masses_copy {
            if mass.applies_gravity{
                grav_masses.push(mass);
            }
        }
        for point in &self.point_masses {
            total_energy_j += point.calculate_total_energy_j(grav_masses.clone());
        }
        
        total_energy_j
    }
    pub fn calculate_system_momentum(&self) -> Vec3 {
        let mut momentum = Vec3::default();
        for point in &self.point_masses {
            momentum += point.calculate_momentum();

        }
        momentum
        
    }
    pub fn draw_historical_positions(&self, window : &mut Window) {
        for point in &self.point_masses{
            point.draw_historical_positions(window);
        }
        for rope in &self.ropes {
            rope.draw_historical_positions(window);
        }
    
    }


}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn conserve_energy_earth_moon() {
        let mut sim = Simulation::default();
        sim.point_masses.push(PointMass::new_earth(None));
        sim.point_masses.push(PointMass::new_moon(None));

        let simulation_jiffy_s = SimulationFloat::new(5_f32);
        dbg!(sim.calculate_total_energy_j());
        let initial_energy = sim.calculate_total_energy_j();
        for interval in 1..5 {
           let criteria = move |simulation_elapsed_time_s : SimulationFloat| {
                    simulation_elapsed_time_s < SimulationFloat::new(589680_f64*(interval as f64)) 
                };
            sim.run_till(&simulation_jiffy_s, criteria.clone());
            dbg!(sim.calculate_total_energy_j());
            dbg!(sim.calculate_system_momentum());
            let mut drift = sim.calculate_total_energy_j() - initial_energy.clone();
            drift.abs_mut();
            assert!(drift < SimulationFloat::new(100_f64));
            //dbg!(sim.point_masses.clone());

        }
    }
    #[test]
    fn conserve_momentum () {
        let mut sim = Simulation::default();
        sim.point_masses.push(PointMass::new_earth(None));
        sim.point_masses.push(PointMass::new_moon_incline(None));

        let simulation_jiffy_s = SimulationFloat::new(5_f32);
        let initial_momentum = sim.calculate_system_momentum();
        for interval in 1..5 {
           let criteria = move |simulation_elapsed_time_s : SimulationFloat| {
                    simulation_elapsed_time_s < SimulationFloat::new(589680_f64*(interval as f64)) 
                };
            sim.run_till(&simulation_jiffy_s, criteria.clone());
            let drift = sim.calculate_system_momentum() - initial_momentum.clone();
            assert!(drift.modulus() < SimulationFloat::new(1e-10_f64));
        }


    }
    #[test]
    fn moon_apogee() {
        
        let mut sim = Simulation::default();
        sim.point_masses.push(PointMass::new_earth(None));
        sim.point_masses.push(PointMass::new_moon_incline(None));

        let simulation_jiffy_s = SimulationFloat::new(5_f32);
        let mut max_lunar_distance = SimulationFloat::new(0_f32);
        let mut min_velocity = Vec3::default();
        let mut min_speed = SimulationFloat::default();

        let steps = 500;
        for step in 1..steps {
            dbg!(step);
            let oribal_period_s = SimulationFloat::new(2358720_f32);
            let criteria = move |simulation_elapsed_time_s : SimulationFloat| {
                    simulation_elapsed_time_s < (oribal_period_s.clone() * SimulationFloat::new((step as f64)/(steps as f64))) // one week of simulation time
                };
            sim.run_till(&simulation_jiffy_s, criteria.clone());


            let lunar_distance = (sim.point_masses[1].pos.clone() - sim.point_masses[0].pos.clone()).modulus();
            if lunar_distance > max_lunar_distance {
                max_lunar_distance = lunar_distance;
                min_velocity = sim.point_masses[1].vel.clone();
                min_speed = min_velocity.modulus();
            }

        }
        dbg!(max_lunar_distance);
        dbg!(min_velocity);
        dbg!(min_speed);
    }

    
}
