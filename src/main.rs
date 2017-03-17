extern crate rand;

#[macro_use]
extern crate glium;
extern crate image;
extern crate time;

use std::io::Cursor;
use glium::{DisplayBuild, Surface};

use glium::glutin;
use std::cmp::Ordering; // needed for wall intersection struct

const MAZE_SIZE: i32 = 5;
const PLAYER_HEIGHT: i32 = 32;
const CELL_SIZE: i32 = 64;
const PROJPLANE_HEIGHT: i32 = 200;
const PROJPLANE_WIDTH: i32 = 320;
//const PLAYER_FOV: f32 = 1.0472; // 60 degrees in radians.

const ANGLE_60: i32 = PROJPLANE_WIDTH;
const ANGLE_30: i32 = ANGLE_60 / 2;
const ANGLE_90: i32 = ANGLE_30 * 3;
const ANGLE_180: i32 = ANGLE_90 * 2;
const ANGLE_270: i32 = ANGLE_90 * 3;
const ANGLE_360: i32 = ANGLE_60 * 6;
const ANGLE_0: i32 = 0;

const PLAYER_FOV: i32 = ANGLE_60;

#[derive(Clone, Copy)]
struct WallIntersection {
    distance_comp: i32,
    distance: f32,
    height: f32,
    //ceiling_height: f32,
    texture_offset: u32,
    wall_intersection_x: f32,
    wall_intersection_y: f32,
}

impl Ord for WallIntersection {
    fn cmp(&self, other: &WallIntersection) -> Ordering {
        self.distance_comp.cmp(&other.distance_comp)
    }
}

impl PartialOrd for WallIntersection {
    fn partial_cmp(&self, other: &WallIntersection) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for WallIntersection {
    fn eq(&self, other: &WallIntersection) -> bool {
        self.distance_comp == other.distance_comp
    }
}

impl Eq for WallIntersection {}

struct TrigTables {
    sin_table: Vec<f32>,
    i_sin_table: Vec<f32>,
    cos_table: Vec<f32>,
    i_cos_table: Vec<f32>,
    tan_table: Vec<f32>,
    i_tan_table: Vec<f32>,
    fish_table: Vec<f32>,
    x_step_table: Vec<f32>,
    y_step_table: Vec<f32>,
}

fn arc_to_rad(arc: f32) -> f32 {
    return (arc * std::f32::consts::PI) / (ANGLE_180 as f32);
}

fn create_tables() -> TrigTables {
    let mut tables = TrigTables {
        sin_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        i_sin_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        cos_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        i_cos_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        tan_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        i_tan_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        fish_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        x_step_table: Vec::with_capacity(ANGLE_360 as usize + 1),
        y_step_table: Vec::with_capacity(ANGLE_360 as usize + 1),
    };

    for i in 0usize..ANGLE_360 as usize {
        // Without the + 0.00001 there will be visual artifacts with rays that intersect walls completely perpendicular.
        let rad = arc_to_rad(i as f32) + 0.0001f32;
        tables.sin_table.push(rad.sin());
        tables.i_sin_table.push(1.0/tables.sin_table[i]);
        tables.cos_table.push(rad.cos());
        tables.i_cos_table.push(1.0/tables.cos_table[i]);
        tables.tan_table.push(rad.tan());
        tables.i_tan_table.push(1.0/tables.tan_table[i]);

        // facing left
        if i >= ANGLE_90 as usize && i < ANGLE_270 as usize
        {
            tables.x_step_table.push(CELL_SIZE as f32/tables.tan_table[i]);
            if tables.x_step_table[i] > 0f32 {
                tables.x_step_table[i] = -tables.x_step_table[i];
            }
        }
        else {
            tables.x_step_table.push(CELL_SIZE as f32/tables.tan_table[i]);
            if tables.x_step_table[i] < 0f32 {
                tables.x_step_table[i] = -tables.x_step_table[i];
            }
        }
        // facing down
        if i >= ANGLE_0 as usize && i < ANGLE_180 as usize {
            tables.y_step_table.push(CELL_SIZE as f32 * tables.tan_table[i]);
            if tables.y_step_table[i] < 0f32 {
                tables.y_step_table[i] = -tables.y_step_table[i];
            }
        }
        else {
            tables.y_step_table.push(CELL_SIZE as f32 * tables.tan_table[i]);
            if tables.y_step_table[i] > 0f32 {
                tables.y_step_table[i] = -tables.y_step_table[i];
            }
        }
    }

    for i in (-ANGLE_30)..ANGLE_30 {
        let rad = arc_to_rad(i as f32);

        tables.fish_table.push(1.0/rad.cos());
    }

    tables
}

// TODO: pack similar data into structs so it's easier to pass this stuff around...
fn draw_floor_column(pixel_x : i32,
                    pixel_y : i32,
                    bottom_y: i32,
                    wall_distance : f32,
                    wall_intersection_x : f32,
                    wall_intersection_y : f32,
                    cell_height : f32,
                    distance_to_projplane : f32,
                    projection_plane_center_y : i32,
                    player_x : f32,
                    player_y : f32,
                    image_width : u32,
                    image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                    dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) {
    // Floor Casting
    let mut current_pixel_y = std::cmp::min(PROJPLANE_HEIGHT-1, pixel_y);

    let floor_distance_from_eyes = PLAYER_HEIGHT as f32 - cell_height;

    while current_pixel_y >= bottom_y && current_pixel_y > 0
    {
        let distance_from_center_to_proj_pixel = projection_plane_center_y - current_pixel_y;
        let staight_distance_from_player_to_floor_intersect = (floor_distance_from_eyes / distance_from_center_to_proj_pixel as f32) * distance_to_projplane as f32;
        // Take the relative angle between the player_angle and the cast ray
        let distance_from_player_to_floor_intersect = staight_distance_from_player_to_floor_intersect;// * beta; //TODO: get to work with proper distance.

        let weight = distance_from_player_to_floor_intersect / wall_distance;

        let floor_intersection_x = weight * wall_intersection_x + (1.0f32 - weight) * player_x;
        let floor_intersection_y = weight * wall_intersection_y + (1.0f32 - weight) * player_y;

        let floor_texture_offset_x = floor_intersection_x as u32 % CELL_SIZE as u32;
        let floor_texture_offset_y = floor_intersection_y as u32 % CELL_SIZE as u32;
        let tex_x = (floor_texture_offset_x as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;
        let tex_y = (floor_texture_offset_y as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

        dest_buffer[current_pixel_y as usize][pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        current_pixel_y -= 1;
    }
}

fn draw_ceiling_column(pixel_x : i32,
                        pixel_y : i32,
                        wall_distance : f32,
                        wall_intersection_x : f32,
                        wall_intersection_y : f32,
                        cell_height : f32,
                        distance_to_projplane : f32,
                        projection_plane_center_y : i32,
                        player_x : f32,
                        player_y : f32,
                        image_width : u32,
                        image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                        dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) {
    // Floor Casting
    let mut current_pixel_y = std::cmp::max(0, pixel_y);

    let floor_distance_from_eyes = cell_height - PLAYER_HEIGHT as f32;

    while current_pixel_y < PROJPLANE_HEIGHT
    {
        let distance_from_center_to_proj_pixel = current_pixel_y - projection_plane_center_y;
        let staight_distance_from_player_to_floor_intersect = (floor_distance_from_eyes / distance_from_center_to_proj_pixel as f32) * distance_to_projplane as f32;
        // Take the relative angle between the player_angle and the cast ray
        let distance_from_player_to_floor_intersect = staight_distance_from_player_to_floor_intersect;// * beta.cos(); TODO: get to work with proper distance.

        let weight = distance_from_player_to_floor_intersect / wall_distance;

        let floor_intersection_x = weight * wall_intersection_x + (1.0f32 - weight) * player_x;
        let floor_intersection_y = weight * wall_intersection_y + (1.0f32 - weight) * player_y;

        let floor_texture_offset_x = floor_intersection_x as u32 % CELL_SIZE as u32;
        let floor_texture_offset_y = floor_intersection_y as u32 % CELL_SIZE as u32;
        let tex_x = (floor_texture_offset_x as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;
        let tex_y = (floor_texture_offset_y as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

        dest_buffer[current_pixel_y as usize][pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        // Ceiling is the same, just mirrored.
        //let ceiling_y = PROJPLANE_HEIGHT as usize - current_pixel_y as usize - 1usize; // Subtract one because it's zero indexed.
        //dest_buffer[ceiling_y][pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        current_pixel_y += 1;
    }
}

fn draw_wall_column(starting_y_coord : i32,
                    column_height : i32,
                    column_bottom : i32,
                    tex_x : usize,
                    image_height : u32,
                    column : i32,
                    image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                    dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) -> i32 {

    let mut proj_y = starting_y_coord;

    if column_height > 0 {
        while proj_y < column_height + column_bottom {
            let tex_y = ((proj_y - column_bottom) as f32 / (column_height as f32 / image_height as f32)) as u32;
            if proj_y >= 0 && proj_y < PROJPLANE_HEIGHT {
                dest_buffer[proj_y as usize][column as usize] = image_buffer[tex_y as usize % image_height as usize][tex_x];
            }
            proj_y += 1;
        }
    }

    proj_y
}

fn main() {
    // building the display, ie. the main object
    let display = glutin::WindowBuilder::new()
        .with_vsync()
        .build_glium()
        .unwrap();

    // building a texture with "OpenGL" drawn on it
    let image = image::load(Cursor::new(&include_bytes!("floor_sample.png")[..]),
                            image::PNG).unwrap().to_rgba();
    let image_dimensions = image.dimensions();
    let (image_width, image_height) = image_dimensions;
    let mut image_raw = image.into_raw().clone();

    // Reversing the image data arranges the pixels correctly but reverses the rgba values of each so they must be reversed again when stored in the image buffer.
    image_raw.reverse();

    let mut image_buffer:Vec< Vec<(u8,u8,u8)> > = Vec::with_capacity(image_height as usize);

    for row_index in 0..image_height {
        let mut row = Vec::with_capacity(image_width as usize);

        for column_index in 0..image_width {

            let image_pixel = row_index as usize * (image_width * 4) as usize + (column_index * 4) as usize;
            //let image_pixel_a = image_raw[image_pixel];
            let image_pixel_b = image_raw[image_pixel+1];
            let image_pixel_g = image_raw[image_pixel+2];
            let image_pixel_r = image_raw[image_pixel+3];

            row.push((image_pixel_r,image_pixel_g, image_pixel_b));
        }

        image_buffer.push(row);
    }

    // building an empty texture
    let dest_texture = glium::Texture2d::empty_with_format(&display,
                                               glium::texture::UncompressedFloatFormat::U8U8U8U8,
                                               glium::texture::MipmapsOption::NoMipmap,
                                               PROJPLANE_WIDTH as u32, PROJPLANE_HEIGHT as u32).unwrap();

    dest_texture.as_surface().clear_color(0.0, 0.0, 0.0, 1.0);

    let mut dest_buffer:Vec< Vec<(u8,u8,u8)> > = Vec::with_capacity(PROJPLANE_HEIGHT as usize);

    for _ in 0..PROJPLANE_HEIGHT {
        let mut row = Vec::with_capacity(PROJPLANE_WIDTH as usize);

        // Init this pixel buffer to greeeeen!!
        for _ in 0..PROJPLANE_WIDTH {
            row.push((0u8,255u8,0u8));
        }

        dest_buffer.push(row);
    }

    // RAYCAST STUFFF!! PEW PEW
    let trig_tables = create_tables();

    let mut player_x = 160f32;
    let mut player_y = 160f32;
    //let mut player_movement_forward = 0f32;
    //let player_turn_speed = (360f32).to_radians();
    let player_move_speed = 100f32;
    let mut player_angle = ANGLE_0 as f32;
    let mut player_arc = player_angle as i32;


    let degree_between_columns = PLAYER_FOV / PROJPLANE_WIDTH;

    let distance_to_projplane = 277f32;//(PROJPLANE_WIDTH as f32 / 2.0) / (PLAYER_FOV as f32 / 2.0).tan();

    // Contains the height of each floor tile.
    let maze = [16f32, 16f32, 16f32, 16f32, 16f32,// 0 - 16f32
                16f32, 0f32, 0f32, 0f32, 16f32, // 32f32 - 128
                16f32, 0f32, 64f32, 0f32, 16f32, // 128 - 196
                16f32, 0f32, 0f32, 0f32, 16f32, // 196 - 256
                16f32, 16f32, 16f32, 16f32, 16f32];

    /*let maze_ceiling = [42, 42, 42, 42, 42,// 0 - 64
                        42, 42, 42, 42, 42, // 64 - 128
                        42, 42, 42, 42, 42, // 128 - 196
                        42, 42, 42, 42, 42, // 196 - 256
                        42, 42, 42, 42, 42];*/


    //let mut last_time = time::precise_time_s() as f32;
    let target_frame_time = 1f32/60f32;
    let delta_time = target_frame_time;

    // Input
    let mut forward_button = false;
    let mut backward_buttom = false;
    let mut turn_left = false;
    let mut turn_right = false;
    let mut cell_draw_limit = 1;

    let mut projection_plane_y_center = PROJPLANE_HEIGHT/2;

    'main: loop {
        { // Raycast rendering
            let player_cell_x = player_x as i32/ CELL_SIZE;
            let player_cell_y = player_y as i32/ CELL_SIZE;

            // Start with the left side of the players view.
            let mut cast_arc = player_arc - (PLAYER_FOV / 2);

            // Wrap
            if cast_arc < 0 {
                cast_arc = ANGLE_360 + cast_arc;
            }

            for column in 0..PROJPLANE_WIDTH {
                //let ray_angle_degrees = ray_anIgle.to_degrees() as i32;
                let mut current_height = maze[(player_cell_y as i32 * MAZE_SIZE + player_cell_x as i32) as usize] as f32;


                // Initially the top of the last wall projected will be the bottom of the projection plane.
                let mut top_of_last_wall = 0i32;

                // horizontal wall intersection data for raycast.
                let mut horz_wall_intersection_y : f32;
                let horz_wall_y_offset : f32;
                let mut horz_wall_intersection_x : f32;

                // vertical wall intersection data for raycast.
                let mut vert_wall_intersection_y : f32;
                let vert_wall_x_offset : f32;
                let mut vert_wall_intersection_x : f32;


                let mut cell_x_offset = 0;
                let mut cell_y_offset = 0;

                // Facing up
                // Reminder coordinate system has y growing downward and x growing to the left
                if cast_arc > ANGLE_0 && cast_arc < ANGLE_180 {
                    // Point A will boarder the cell below
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    horz_wall_y_offset = CELL_SIZE as f32;

                    let x_dist_to_intersection = (horz_wall_intersection_y - player_y) as f32 * trig_tables.i_tan_table[cast_arc as usize];
                    horz_wall_intersection_x = x_dist_to_intersection + player_x as f32;
                } else {
                    //  Point A will boarder the cell above
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32;
                    horz_wall_y_offset = -CELL_SIZE as f32;
                    cell_y_offset = -1;

                    let x_dist_to_intersection = (horz_wall_intersection_y - player_y) as f32 * trig_tables.i_tan_table[cast_arc as usize];
                    horz_wall_intersection_x = x_dist_to_intersection + player_x as f32;
                }

                // Facing right
                if cast_arc < ANGLE_90 || cast_arc > ANGLE_270
                {
                    // Point A will boarder the cell to the right
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    vert_wall_x_offset = CELL_SIZE as f32;

                    let y_dist_to_intersection : f32 = (vert_wall_intersection_x - player_x) as f32 * trig_tables.tan_table[cast_arc as usize];
                    vert_wall_intersection_y = y_dist_to_intersection + player_y as f32;
                }
                else
                {
                    //  Point A will boarder the cell to the left
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE)  * CELL_SIZE) as f32;
                    vert_wall_x_offset = -CELL_SIZE as f32;
                    cell_x_offset = -1;

                    let y_dist_to_intersection : f32 = (vert_wall_intersection_x - player_x) as f32 * trig_tables.tan_table[cast_arc as usize];
                    vert_wall_intersection_y = y_dist_to_intersection + player_y as f32;
                }


                let mut horz_out_of_range = false;
                let mut vert_out_of_range = false;
                let mut wall_intersections: Vec<WallIntersection> = Vec::new();
                while !horz_out_of_range || !vert_out_of_range {
                    if cast_arc != ANGLE_180 && cast_arc != ANGLE_0 {
                        let horz_wall_x_offset = trig_tables.x_step_table[cast_arc as usize];
                        loop {
                            // Will be used to index into the map array to check for walls.
                            let cell_y : i32 = horz_wall_intersection_y as i32 / CELL_SIZE + cell_y_offset;
                            let cell_x : i32 = horz_wall_intersection_x as i32 / CELL_SIZE;
                            // Check the horz_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                            // This will then display a horz wall at the 0 index outside the map.
                            if cell_x < 0i32 || cell_y < 0i32 || cell_x >= MAZE_SIZE || cell_y >= MAZE_SIZE {
                                horz_out_of_range = true;
                                break;
                            } else {
                                let hor_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize];
                                //let ceiling_height = maze_ceiling[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let dist_to_horz_intersection = (horz_wall_intersection_x - player_x) * trig_tables.i_cos_table[cast_arc as usize];

                                let wall_intersection = WallIntersection {
                                    distance_comp: (dist_to_horz_intersection * 1000f32) as i32,
                                    distance: dist_to_horz_intersection,
                                    height: hor_cell_height,
                                    //ceiling_height: ceiling_height,
                                    texture_offset: horz_wall_intersection_x as u32 % CELL_SIZE as u32,
                                    wall_intersection_x: horz_wall_intersection_x,
                                    wall_intersection_y: horz_wall_intersection_y,
                                };

                                wall_intersections.push(wall_intersection);
                            }

                            // Extend the ray further, checking what will be the next intersection.
                            horz_wall_intersection_y += horz_wall_y_offset;
                            horz_wall_intersection_x += horz_wall_x_offset;
                        }
                    }
                    else
                    {
                        horz_out_of_range = true;
                    }

                    if cast_arc != ANGLE_90 && cast_arc != ANGLE_270 {
                        let vert_wall_y_offset = trig_tables.y_step_table[cast_arc as usize];
                        loop {
                            let cell_y = vert_wall_intersection_y as i32 / CELL_SIZE;
                            let cell_x = vert_wall_intersection_x as i32 / CELL_SIZE + cell_x_offset;
                            // Check the vert_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                            // This will then display a vert wall at the 0 index outside the map.
                            if cell_x < 0i32 || cell_y < 0i32 || cell_y >= MAZE_SIZE || cell_x >= MAZE_SIZE {
                                vert_out_of_range = true;
                                break;
                            } else {
                                let vert_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize];
                                //let ceiling_height = maze_ceiling[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let dist_to_vert_intersection = (vert_wall_intersection_y - player_y) * trig_tables.i_sin_table[cast_arc as usize];

                                let wall_intersection = WallIntersection {
                                    distance_comp: (dist_to_vert_intersection * 1000f32) as i32,
                                    distance: dist_to_vert_intersection,
                                    height: vert_cell_height,
                                    //ceiling_height: ceiling_height,
                                    texture_offset: vert_wall_intersection_y as u32 % CELL_SIZE as u32,
                                    wall_intersection_x: vert_wall_intersection_x,
                                    wall_intersection_y: vert_wall_intersection_y,
                                };

                                wall_intersections.push(wall_intersection);
                            }

                            // Extend the ray further, checking what will be the next intersection.
                            vert_wall_intersection_y += vert_wall_y_offset;
                            vert_wall_intersection_x += vert_wall_x_offset;
                        }
                    }
                    else
                    {
                        vert_out_of_range = true;
                    }
                }

                // Sort all the wall intersections such that the closer walls will be drawn first.
                wall_intersections.sort();

                // Draw walls and floors.
                let mut draw_count = 0;
                for intersection in &wall_intersections {
                    /*if draw_count >= cell_draw_limit {
                        break;
                    }*/

                    // The height - the change in height will be the portion of the wall hidden (floor/top of the last wall will cover it)
                    let last_height = current_height;

                    // Continous slices can be drawn all at once.
                    if last_height == intersection.height {
                        continue;
                    }

                    current_height = intersection.height;

                    let sunken = last_height > current_height;

                    // Correct "fishbowl effect". Fish table angle is the angle of the cast ray, relative to the player's viewing angle.
                    let distance = intersection.distance / trig_tables.fish_table[column as usize];

                    let last_projected_slice_height = (last_height as f32 / distance * distance_to_projplane) as i32;
                    let projected_slice_height = (current_height / distance * distance_to_projplane) as i32;

                    // TODO: write the formula in long form in comments here.
                    // TODO: The further from 32 (Player's eye level), there appears to be some rounding woes that cause the wall to render awkwardly (jumpy?)
                    let divisor = distance / distance_to_projplane;
                    let projected_bottom_coord = if divisor > std::f32::EPSILON {
                        projection_plane_y_center - (PLAYER_HEIGHT as f32 / divisor) as i32
                    } else {
                        0
                    };

                    //let projected_bottom_coord = projection_plane_y_center - (projected_slice_height as f32 * 0.5f32) as i32;

                    // Start rendering from height of the last drawn wall.
                    let starting_y_coord = projected_bottom_coord + last_projected_slice_height;
                    // Subtract the starting_y_coord by 1 since that starting y is for the wall, while the floor will be drawn starting from one below it.
                    let (projplane_pixel_x, projplane_pixel_y) = (column, starting_y_coord-1);

                    // TODO: have different textures, super confusing right now.
                    draw_floor_column(projplane_pixel_x,
                                    projplane_pixel_y,
                                    top_of_last_wall,
                                    distance,
                                    intersection.wall_intersection_x,
                                    intersection.wall_intersection_y,
                                    last_height,
                                    distance_to_projplane,
                                    projection_plane_y_center,
                                    player_x,
                                    player_y,
                                    image_width,
                                    &image_buffer,
                                    &mut dest_buffer);

                    // No need to draw the wall when it is sunken (as it will not be visible).
                    //if sunken {
                        //top_of_last_wall = starting_y_coord;
                    //} else {
                        let tex_x = (intersection.texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

                        // Calculate where the rendering of a wall should start.
                        let wall_bottom_coord = if starting_y_coord > top_of_last_wall { starting_y_coord } else { top_of_last_wall };
                        top_of_last_wall = draw_wall_column(wall_bottom_coord, projected_slice_height, projected_bottom_coord, tex_x, image_height, column, &image_buffer, &mut dest_buffer);
                    //}

                    /*draw_ceiling_column(projplane_pixel_x,
                                        top_of_last_wall,
                                        distance,
                                        intersection.wall_intersection_x,
                                        intersection.wall_intersection_y,
                                        current_height,
                                        distance_to_projplane,
                                        projection_plane_y_center,
                                        player_x,
                                        player_y,
                                        image_width,
                                        &image_buffer,
                                        &mut dest_buffer);*/

                    draw_count += 1;
                }

                cast_arc += degree_between_columns;

                if cast_arc >= ANGLE_360 {
                    cast_arc -= ANGLE_360;
                }
            }
        }

        // drawing a frame
        let target = display.draw();
        let the_texture = glium::Texture2d::new(&display, dest_buffer.clone()).unwrap();
        the_texture.as_surface().fill(&target, glium::uniforms::MagnifySamplerFilter::Nearest);
        target.finish().unwrap();

        // Clear Buffer
        dest_texture.as_surface().clear_color(0.0, 0.0, 0.0, 1.0);
        for row in 0..PROJPLANE_HEIGHT {
            for cell in 0..PROJPLANE_WIDTH {
                dest_buffer[row as usize][cell as usize] = (0u8, 255u8, 0u8);
            }
        }

        // TODO: Store the buttons which have been pressed in a map or array.
        for ev in display.poll_events() {
            match ev {
                glium::glutin::Event::Closed => return,
                glium::glutin::Event::KeyboardInput(element_state, _, virtual_key_code) => {
                    match element_state {
                        glium::glutin::ElementState::Pressed => {
                            let key = virtual_key_code.unwrap();
                            match key {
                                glium::glutin::VirtualKeyCode::W => {
                                    forward_button = true;
                                },
                                glium::glutin::VirtualKeyCode::A => {
                                    turn_left = true;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    backward_buttom = true;
                                },
                                glium::glutin::VirtualKeyCode::D => {
                                    turn_right = true;
                                },

                                glium::glutin::VirtualKeyCode::Escape => {
                                    break 'main;
                                }
                                _ => ()
                            }
                        },
                        glium::glutin::ElementState::Released => {
                            let key = virtual_key_code.unwrap();
                            match key {
                                glium::glutin::VirtualKeyCode::W => {
                                    forward_button = false;
                                },
                                glium::glutin::VirtualKeyCode::A => {
                                    turn_left = false;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    backward_buttom = false;
                                },
                                glium::glutin::VirtualKeyCode::D => {
                                    turn_right = false;
                                },
                                glium::glutin::VirtualKeyCode::E => {
                                    cell_draw_limit += 1;
                                },
                                glium::glutin::VirtualKeyCode::Q => {
                                    cell_draw_limit -= 1;
                                },
                                glium::glutin::VirtualKeyCode::R => {
                                    projection_plane_y_center += 5;
                                },
                                glium::glutin::VirtualKeyCode::F => {
                                    projection_plane_y_center -= 5;
                                },
                                _ => ()
                            }
                        }
                    }
                },
                _ => ()
            }
        }

        if turn_left {
            player_angle -= ANGLE_180 as f32 * delta_time;
            if player_angle < ANGLE_0 as f32 {
                player_angle += ANGLE_360 as f32;
            }
        }

        if turn_right {
            player_angle += ANGLE_180 as f32 * delta_time;
            if player_angle >= ANGLE_360 as f32 {
                player_angle -= ANGLE_360 as f32;
            }
        }

        player_arc = player_angle as i32;

        let player_x_dir = trig_tables.cos_table[player_arc as usize];
        let player_y_dir = trig_tables.sin_table[player_arc as usize];

        if forward_button {
            player_x += player_move_speed * delta_time * player_x_dir;
            player_y += player_move_speed * delta_time * player_y_dir;
        }

        if backward_buttom {
            player_x -= player_move_speed * delta_time * player_x_dir;
            player_y -= player_move_speed * delta_time * player_y_dir;
        }

        std::thread::sleep(std::time::Duration::from_millis(20));
    }
}
