extern crate rand;

#[macro_use]
extern crate glium;
extern crate image;
extern crate time;

use std::io::Cursor;
use glium::{DisplayBuild, Surface};

use glium::glutin;
use std::cmp::Ordering; // needed for wall intersection struct

const MAZE_SIZE : i32 = 5;
const PLAYER_HEIGHT : i32 = 32;
const CELL_SIZE : i32 = 64;
const PROJPLANE_HEIGHT : i32 = 200;
const PROJPLANE_WIDTH : i32 = 320;
const PLAYER_FOV : f32 = 1.0472; // 60 degrees in radians.

#[derive(Clone, Copy, Eq)]
struct WallIntersection {
    distance: i32,
    height: i32,
    texture_offset: u32,
    wall_intersection_x: i32,
    wall_intersection_y: i32
}

impl Ord for WallIntersection {
    fn cmp(&self, other: &WallIntersection) -> Ordering {
        self.distance.cmp(&other.distance)
    }
}

impl PartialOrd for WallIntersection {
    fn partial_cmp(&self, other: &WallIntersection) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for WallIntersection {
    fn eq(&self, other: &WallIntersection) -> bool {
        self.distance == other.distance
    }
}

// TODO: pack similar data into structs so it's easier to pass this stuff around...
fn draw_floor_column(pixel_x : i32,
                    pixel_y : i32,
                    bottom_y: i32,
                    wall_distance : f32,
                    wall_intersection_x : f32,
                    wall_intersection_y : f32,
                    distance_to_projplane : f32,
                    player_x : f32,
                    player_y : f32,
                    image_width : u32,
                    image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                    dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) {
    // Floor Casting
    let mut current_pixel_y = pixel_y;
    let projection_plane_center_y = PROJPLANE_HEIGHT/2;

    current_pixel_y = std::cmp::min(PROJPLANE_HEIGHT-1, pixel_y);

    while current_pixel_y >= bottom_y && current_pixel_y > 0
    {
        let distance_from_center_to_proj_pixel = projection_plane_center_y - current_pixel_y;
        let staight_distance_from_player_to_floor_intersect = (PLAYER_HEIGHT as f32 / distance_from_center_to_proj_pixel as f32) * distance_to_projplane as f32;
        // Take the relative angle between the player_angle and the cast ray
        let distance_from_player_to_floor_intersect = staight_distance_from_player_to_floor_intersect; //* beta.cos(); TODO: get to work with proper distance.

        let weight = distance_from_player_to_floor_intersect / wall_distance;

        let floor_intersection_x = weight * wall_intersection_x + (1.0f32 - weight) * player_x;
        let floor_intersection_y = weight * wall_intersection_y + (1.0f32 - weight) * player_y;

        let floor_texture_offset_x = floor_intersection_x as u32 % CELL_SIZE as u32;
        let floor_texture_offset_y = floor_intersection_y as u32 % CELL_SIZE as u32;
        let tex_x = (floor_texture_offset_x as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;
        let tex_y = (floor_texture_offset_y as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

        dest_buffer[current_pixel_y as usize][pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        //     // Ceiling is the same, just mirrored.
        //     let ceiling_y = PROJPLANE_HEIGHT as usize - projplane_pixel_y as usize - 1usize; // Subtract one because it's zero indexed.
        //     dest_buffer[ceiling_y][projplane_pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        current_pixel_y -= 1;
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
            // TODO: Look into this method of using ints to avoid floats.
            let d = proj_y * 256 - PROJPLANE_HEIGHT * 128 + column_height * 128;  //256 and 128 factors to avoid floats
            let tex_y = ((d * image_height as i32) / column_height) / 256;
            if proj_y >= 0 && proj_y < PROJPLANE_HEIGHT {
                dest_buffer[proj_y as usize][column as usize] = image_buffer[tex_y as usize % image_height as usize][tex_x];
            }
            proj_y += 1;
        }
    }

    proj_y
}

fn draw_wall_and_floor_columns(player_x: f32,
                            player_y: f32,
                            wall_intersection_x: i32,
                            wall_intersection_y: i32,
                            cell_height: f32,
                            distance: f32,
                            distance_to_projplane: f32,
                            column: i32,
                            top_of_last_wall: i32,
                            wall_texture_offset: u32,
                            image_width: u32,
                            image_height: u32,
                            image_buffer: &Vec< Vec<(u8,u8,u8)> >,
                            mut dest_buffer: &mut Vec< Vec<(u8,u8,u8)> >) -> i32 {

    let projected_slice_height = (cell_height / distance * distance_to_projplane) as i32;
    let projected_bottom_coord = PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / (distance / distance_to_projplane)) as i32;

    // Floor Casting TODO: bring back once wall visiual glitches are sorted.
    let (projplane_pixel_x, projplane_pixel_y) = (column, projected_bottom_coord);

    // TODO: have different textures, super confusing right now.
    // TODO: Readd the ceilling drawing.
    /*draw_floor_column(projplane_pixel_x,
                    projplane_pixel_y,
                    top_of_last_wall,
                    distance,
                    wall_intersection_x as f32,
                    wall_intersection_y as f32,
                    distance_to_projplane,
                    player_x,
                    player_y,
                    image_width,
                    &image_buffer,
                    &mut dest_buffer);*/

    let tex_x = (wall_texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

    let wall_bottom_coord = if projected_bottom_coord > top_of_last_wall { projected_bottom_coord } else { top_of_last_wall };

    draw_wall_column(wall_bottom_coord, projected_slice_height, projected_bottom_coord, tex_x, image_height, column, &image_buffer, &mut dest_buffer);

    projected_bottom_coord + projected_slice_height
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
    let mut player_x = 160f32;
    let mut player_y = 160f32;
    let mut player_movement_forward = 0f32;
    let player_turn_speed = (360f32).to_radians();
    let player_move_speed = 100f32;
    let mut player_angle = (0f32).to_radians();

    let degree_between_columns = PLAYER_FOV / PROJPLANE_WIDTH as f32;

    let distance_to_projplane = (PROJPLANE_WIDTH as f32 / 2.0) / (PLAYER_FOV as f32 / 2.0).tan();

    // Contains the height of each floor tile.
    let maze = [72, 72, 72, 72, 72,// 0 - 64
                72, 0, 40, 0, 72, // 64 - 128
                72, 0, 12, 24, 72, // 128 - 196
                72, 0, 0, 0, 72, // 196 - 256
                72, 72, 72, 72, 72];


    //let mut last_time = time::precise_time_s() as f32;
    let target_frame_time = 1f32/60f32;
    let delta_time = target_frame_time;

    // Input
    let mut forward_button = false;
    let mut backward_buttom = false;
    let mut turn_left = false;
    let mut turn_right = false;

    'main: loop {
        { // Raycast rendering
            let player_cell_x = (player_x / CELL_SIZE as f32) as i32;
            let player_cell_y = (player_y / CELL_SIZE as f32) as i32;

            // Start with the left side of the players view. Without the + 0.00001 there will be visual artifacts with rays that intersect walls completely perpendicular.
            let mut ray_angle : f32 = player_angle - (PLAYER_FOV / 2.0) + 0.00001;

            // Wrap
            if ray_angle < 0f32 {
                ray_angle = (360f32).to_radians() + ray_angle;
            }

            for column in 0..PROJPLANE_WIDTH {
                let ray_angle_degrees = ray_angle.to_degrees() as i32;
                let mut current_height = maze[(player_cell_y * MAZE_SIZE + player_cell_x) as usize];

                // Initially the top of the last wall projected will be the bottom of the projection plane.
                let mut top_of_last_wall = 0i32;

                // horizontal wall intersection data for raycast.
                let mut horz_wall_intersection_y : f32;
                let horz_wall_y_offset : f32;
                let mut horz_wall_x_offset : f32;
                let mut horz_wall_intersection_x : f32;

                // vertical wall intersection data for raycast.
                let mut vert_wall_intersection_y : f32;
                let mut vert_wall_y_offset : f32;
                let vert_wall_x_offset : f32;
                let mut vert_wall_intersection_x : f32;

                horz_wall_x_offset = CELL_SIZE as f32 / ray_angle.tan();
                vert_wall_y_offset = CELL_SIZE as f32 * ray_angle.tan();

                // Facing up
                // Reminder coordinate system has y growing downward and x growing to the left
                if ray_angle >= (0f32).to_radians() && ray_angle < (180f32).to_radians() {
                    if vert_wall_y_offset < 0.0 {
                        vert_wall_y_offset *= -1f32;
                    }
                    // Point A will boarder the cell below
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    horz_wall_y_offset = CELL_SIZE as f32;
                } else {
                    if vert_wall_y_offset > 0.0 {
                        vert_wall_y_offset *= -1f32;
                    }
                    //  Point A will boarder the cell above
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32;
                    horz_wall_y_offset = -CELL_SIZE as f32;
                }

                // TODO: Handle when the angle is zero/360 or 180 degrees (parallel to horizontal walls)
                let mut x_dist_to_intersection : f32 = std::f32::MAX;
                if ray_angle_degrees != 0 && ray_angle_degrees != 180 {
                    x_dist_to_intersection = (horz_wall_intersection_y - player_y) as f32 / ray_angle.tan();
                }

                horz_wall_intersection_x = x_dist_to_intersection + player_x as f32;

                // Facing right
                if ray_angle < (90f32).to_radians() || ray_angle >= (270f32).to_radians()
                {
                    if horz_wall_x_offset < 0f32 {
                        horz_wall_x_offset *= -1f32;
                    }
                    // Point A will boarder the cell to the right
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    vert_wall_x_offset = CELL_SIZE as f32;
                }
                else
                {
                    if horz_wall_x_offset > 0f32 {
                        horz_wall_x_offset *= -1f32;
                    }
                    //  Point A will boarder the cell to the left
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE)  * CELL_SIZE) as f32; // TODO: try moving this -1 maybe it'll help fix the corner issue (def is having an effect).
                    vert_wall_x_offset = -CELL_SIZE as f32;
                }

                let y_dist_to_intersection : f32 = (vert_wall_intersection_x - player_x) as f32 * ray_angle.tan();
                vert_wall_intersection_y = y_dist_to_intersection + player_y as f32;

                // TODO: remove the rounded ray_angle_degrees?
                if !(ray_angle >= (0f32).to_radians() && ray_angle < (180f32).to_radians()) {
                    horz_wall_intersection_y -= 1f32;
                }

                if !(ray_angle < (90f32).to_radians() || ray_angle >= (270f32).to_radians()) {
                    vert_wall_intersection_x -= 1f32;
                }

                let mut horz_out_of_range = false;
                let mut vert_out_of_range = false;
                let mut wall_intersections: Vec<WallIntersection> = Vec::new();
                while !horz_out_of_range || !vert_out_of_range {
                    if ray_angle_degrees != 180 && ray_angle_degrees != 0 {
                        loop {
                            // Will be used to index into the map array to check for walls.
                            let cell_y : i32 = (horz_wall_intersection_y / CELL_SIZE as f32) as i32;
                            let cell_x : i32 = (horz_wall_intersection_x / CELL_SIZE as f32) as i32;
                            // Check the horz_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                            // This will then display a horz wall at the 0 index outside the map.
                            if horz_wall_intersection_x < 0f32 || horz_wall_intersection_y < 0f32 || horz_wall_intersection_x >= (MAZE_SIZE * CELL_SIZE) as f32 || horz_wall_intersection_y >= (MAZE_SIZE * CELL_SIZE) as f32 {
                                horz_out_of_range = true;
                                break;
                            } else /*if maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0*/ {
                                let hor_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let y_dist = (horz_wall_intersection_y - player_y).powf(2.0);
                                let x_dist = (horz_wall_intersection_x - player_x).powf(2.0);
                                let dist_to_horz_intersection = (x_dist + y_dist).sqrt();

                                let wall_intersection = WallIntersection {
                                    distance: dist_to_horz_intersection as i32,
                                    height: hor_cell_height as i32,
                                    texture_offset: horz_wall_intersection_x as u32 % CELL_SIZE as u32,
                                    wall_intersection_x: horz_wall_intersection_x as i32,
                                    wall_intersection_y: horz_wall_intersection_y as i32
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

                    if ray_angle_degrees != 90 && ray_angle_degrees != 270 {
                        loop {
                            let cell_y = (vert_wall_intersection_y / CELL_SIZE as f32) as i32;
                            let cell_x = (vert_wall_intersection_x / CELL_SIZE as f32) as i32;
                            // Check the vert_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                            // This will then display a vert wall at the 0 index outside the map.
                            if vert_wall_intersection_x < 0f32 || vert_wall_intersection_y < 0f32 || cell_y >= MAZE_SIZE || cell_x >= MAZE_SIZE {
                                vert_out_of_range = true;
                                break;
                            } else /*if maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0*/ {
                                let vert_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let y_dist = (vert_wall_intersection_y - player_y).powf(2.0);
                                let x_dist = (vert_wall_intersection_x - player_x).powf(2.0);
                                let dist_to_vert_intersection = (x_dist + y_dist).sqrt();

                                let wall_intersection = WallIntersection {
                                    distance: dist_to_vert_intersection as i32,
                                    height: vert_cell_height as i32,
                                    texture_offset: vert_wall_intersection_y as u32 % CELL_SIZE as u32,
                                    wall_intersection_x: vert_wall_intersection_x as i32,
                                    wall_intersection_y: vert_wall_intersection_y as i32
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
                //let mut draw = true;
                for intersection in &wall_intersections {
                    // TODO: I think i've found a problem with only drawing raised and not sunken walls...
                    // Or maybe I need to detect something about there being no walls??
                    // Maybe if I keep track of all the cells and then draw floors also when there are sunken walls detected (raised platform in center of room)
                    // Raised Wall
                    //if intersection.height > current_height {
                        // The height - the change in height will be the portion of the wall hidden (floor/top of the last wall will cover it)
                        let last_height = current_height;
                        current_height = intersection.height;

                        // Correct "fishbowl effect". Beta angle is the angle of the cast ray, relative to the player's viewing angle.
                        // TODO: This can be done once per column instead of per every wall
                        let beta = ray_angle - player_angle;
                        let distance = intersection.distance as f32 * beta.cos();

                        let last_projected_slice_height = (last_height as f32 / distance * distance_to_projplane) as i32;
                        let projected_slice_height = (intersection.height as f32 / distance * distance_to_projplane) as i32;

                        // TODO: write the formula in long form in comments here.
                        let divisor = distance / distance_to_projplane;
                        let projected_bottom_coord = if divisor > std::f32::EPSILON {
                            PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / divisor) as i32
                        } else {
                            0
                        };

                        // Start rendering from height of the last drawn wall.
                        let starting_y_coord = projected_bottom_coord + last_projected_slice_height;

                        let (projplane_pixel_x, projplane_pixel_y) = (column, starting_y_coord);

                        // TODO: have different textures, super confusing right now.
                        // TODO: Readd the ceilling drawing.
                        /*draw = !draw;
                        if !draw {
                            top_of_last_wall = projected_slice_height + projected_bottom_coord;
                            continue;
                        }*/
                        draw_floor_column(projplane_pixel_x,
                                        projplane_pixel_y,
                                        top_of_last_wall,
                                        distance,
                                        intersection.wall_intersection_x as f32,
                                        intersection.wall_intersection_y as f32,
                                        distance_to_projplane,
                                        player_x,
                                        player_y,
                                        image_width,
                                        &image_buffer,
                                        &mut dest_buffer);

                        let tex_x = (intersection.texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

                        // Calculate where the rendering of a wall should start.
                        let wall_bottom_coord = if starting_y_coord > top_of_last_wall { starting_y_coord } else { top_of_last_wall };

                        top_of_last_wall = draw_wall_column(wall_bottom_coord, projected_slice_height, projected_bottom_coord, tex_x, image_height, column, &image_buffer, &mut dest_buffer);
                    //}
                    /*else if intersection.height < current_height
                    {
                        // The height - the change in height will be the portion of the wall hidden (floor/top of the last wall will cover it)
                        let last_height = current_height;
                        current_height = intersection.height;

                        // Correct "fishbowl effect". Beta angle is the angle of the cast ray, relative to the player's viewing angle.
                        // TODO: This can be done once per column instead of per every wall
                        let beta = ray_angle - player_angle;
                        let distance = intersection.distance as f32 * beta.cos();

                        let last_projected_slice_height = (last_height as f32 / distance * distance_to_projplane) as i32;
                        let projected_bottom_coord = PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / (distance / distance_to_projplane)) as i32;

                        // Start rendering from height of the last drawn wall.
                        let starting_y_coord = projected_bottom_coord + last_projected_slice_height;

                        let (projplane_pixel_x, projplane_pixel_y) = (column, starting_y_coord);

                        // TODO: have different textures, super confusing right now.
                        // TODO: Readd the ceilling drawing.
                        draw_floor_column(projplane_pixel_x,
                                        projplane_pixel_y,
                                        top_of_last_wall,
                                        distance,
                                        intersection.wall_intersection_x as f32,
                                        intersection.wall_intersection_y as f32,
                                        distance_to_projplane,
                                        player_x,
                                        player_y,
                                        image_width,
                                        &image_buffer,
                                        &mut dest_buffer);
                    }*/
                }

                ray_angle += degree_between_columns;

                if ray_angle as i32 >= (360f32).to_radians() as i32 {
                    ray_angle = ray_angle - (360f32).to_radians();
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
                                    backward_buttom = true;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    turn_left = true;
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
                                    backward_buttom = false;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    turn_left = false;
                                },
                                glium::glutin::VirtualKeyCode::D => {
                                    turn_right = false;
                                },
                                _ => ()
                            }
                        }
                    }
                },
                _ => ()
            }
        }

        if forward_button {
            player_movement_forward = player_move_speed;
        }
        if backward_buttom {
            player_angle -= player_turn_speed * delta_time;
        }
        if turn_left {
            player_movement_forward = -player_move_speed;
        }
        if turn_right {
            player_angle += player_turn_speed * delta_time;
        }

        player_x += player_movement_forward * delta_time * player_angle.cos();
        player_y += player_movement_forward * delta_time * player_angle.sin();
        player_movement_forward = 0f32;

        // Wrap
        if player_angle < 0f32 {
            player_angle = (360f32).to_radians() + player_angle;
        } else if player_angle >= (360f32).to_radians() {
            player_angle = player_angle - (360f32).to_radians();
        }

        std::thread::sleep(std::time::Duration::from_millis(16));
    }
}
