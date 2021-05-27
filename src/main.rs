use clap::{App, load_yaml};
use libnotcurses_sys::*;
use bio::io::fasta::Reader;

extern crate ndarray;

use ndarray::prelude::*;
use std::cmp::max;

use std::collections::BTreeMap;

struct PolyRecord
{
	genes: Vec<String>,
	records: Array2<u8>,
	max_name_length: u32,
	max_sequence_length: u32,
	top_id: usize,
}

impl PolyRecord
{
	fn new(file_name: &str, primary_sequence: &str) -> Option<PolyRecord>
	{
		let fasta_file = Reader::from_file(file_name);

		let fasta_file = match fasta_file
		{
			Ok(db) => db,
			Err(_error) => return None,
		};

		let mut genes = Vec::new();
		let mut records = Array2::<u8>::zeros((0, 0));

		let mut top_id: usize = 0;

		let mut max_sequence_length: u32 = 0;
		let mut max_name_length: u32 = 0;

		for (index, record) in fasta_file.records().enumerate()
		{
			let record = record.unwrap();

			if records.ncols() == 0
			{
				records = Array2::<u8>::zeros((0, record.seq().len()));
			}

			if record.id() == primary_sequence || (primary_sequence == "" && index == 0)
			{
				top_id = index;
			}

			records.push_row(ArrayView::from(&record.seq()));
			genes.push(record.id().to_string());

			max_name_length = max(max_name_length, record.id().len() as u32);
			max_sequence_length = max(max_sequence_length, record.seq().len() as u32);
		}

		Some(PolyRecord
		{
			genes,
			records,
			max_name_length,
			max_sequence_length,
			top_id,
		})
	}

	fn draw_id(&self, plane: &mut NcPlane, id_position: usize)
	{
		let id_name = &self.genes[id_position];

		let difference: u32 = self.max_name_length - id_name.len() as u32;
		if difference > 0
		{
			plane.cursor_move_cols(difference as i32);
		}

		plane.putstrln(&id_name);
	}

	fn draw_sequence(&self, plane: &mut NcPlane, sequence: &Array1<u8>, position: usize)
	{
		let mut tmp = [0; 4];

		let top_record = &self.records.row(self.top_id as usize);

		for (index, letter) in sequence.iter().enumerate()
		{
			if index == 0 || letter != sequence.iter().nth(index - 1).unwrap()
			{
				match letter
				{
					b'A' => { plane.set_bg_rgb8(216, 27, 96); },
					b'T' => { plane.set_bg_rgb8(30, 136, 229); },
					b'G' => { plane.set_bg_rgb8(255, 193, 7); },
					b'C' => { plane.set_bg_rgb8(0, 77, 64); },
					b'-' => { plane.set_bg_rgb8(254, 97, 0); },
					_ => {},
				};
			}

			let mut modified_letter = *letter as char;

			if self.top_id != position && letter != top_record.iter().nth(index).unwrap()
			{
				modified_letter.make_ascii_lowercase();
			}

			plane.putstr(modified_letter.encode_utf8(&mut tmp));
		}

		plane.cursor_move_yx(plane.cursor_y() + 1, 0);
	}

	fn draw(&self) -> NcResult<()>
	{
		let mut tmp = [0; 4];

		let mut nc = Notcurses::with_flags
		(
			NCOPTION_SUPPRESS_BANNERS | NCOPTION_NO_WINCH_SIGHANDLER | NCOPTION_NO_QUIT_SIGHANDLERS,
		)?;

		let (t_rows, _t_cols) = nc.term_dim_yx();
		let mut stdplane = nc.stdplane();

		let mut top_sequences_plane = NcPlane::new_bound(&mut stdplane, 0, self.max_name_length as i32 + 2, 3, self.max_sequence_length)?;
		let mut top_left_plane = NcPlane::new_bound(&mut stdplane, 0, 0, 3, self.max_name_length + 2)?;

		top_left_plane.set_base(" ", 0, NcChannelPair::with_default())?;
		top_sequences_plane.set_base(" ", 0, NcChannelPair::with_default())?;

		let mut temp_position: i32 = 19;
		while temp_position < self.max_sequence_length as i32
		{
			top_sequences_plane.cursor_move_yx(0, temp_position as u32);
			top_sequences_plane.putstr(&(temp_position + 1).to_string());
			top_sequences_plane.cursor_move_yx(1, temp_position as u32);
			top_sequences_plane.putstr("|");

			temp_position += 20;
		}

		top_sequences_plane.cursor_move_yx(2, 0);
		top_sequences_plane.set_fg_rgb8(0, 0, 0);

		top_left_plane.cursor_move_rows(2)?;

		let mut sequences_plane = NcPlane::new_bound(&mut stdplane, 3, self.max_name_length as i32 + 2, self.records.nrows() as u32 - 1, self.max_sequence_length)?;
		let mut left_plane = NcPlane::new_bound(&mut stdplane, 3, 0, self.records.nrows() as u32 - 1, self.max_name_length + 2)?;
		let histogram_plane = NcPlane::new_bound(&mut stdplane, t_rows as i32 - 7, self.max_name_length as i32 + 2, 20, self.max_sequence_length)?;
		let conserved_plane = NcPlane::new_bound(&mut stdplane, t_rows as i32 - 1, self.max_name_length as i32 + 2, 20, self.max_sequence_length)?;
		let histogram_status_plane = NcPlane::new_bound(&mut stdplane, t_rows as i32 - 7, 0, 20, self.max_name_length + 2)?;

		histogram_status_plane.set_base(" ", 0, NcChannelPair::with_default())?;
		histogram_status_plane.cursor_move_rows(2 as i32)?;
		histogram_plane.set_bg_rgb8(0, 77, 64);

		let difference = self.max_name_length - 10;
		if difference > 0
		{
			histogram_status_plane.cursor_move_cols(difference as i32)?;
		}

		histogram_status_plane.putstr("Conserved")?;
		histogram_status_plane.putstrln(" 1")?;

		let difference = self.max_name_length - 8;
		if difference > 0
		{
			histogram_status_plane.cursor_move_cols(difference as i32)?;
		}

		histogram_status_plane.putstr("Details")?;
		histogram_status_plane.cursor_move_rows(2 as i32)?;
		histogram_status_plane.putstrln(" 0")?;

		left_plane.set_base(" ", 0, NcChannelPair::with_default())?;
		sequences_plane.set_base(".", 0, NcChannelPair::with_rgb(0x000000, 0xFF00))?;

		let mut index = 0;
		for column in self.records.columns()
		{
			let mut counts = BTreeMap::new();
			for letter in column.iter()
			{
			    *counts.entry(letter).or_insert(0) += 1;
			}

			let max_letter = counts.into_iter().max_by_key(|&(_, count)| count).unwrap().0;
			let frequency = column.iter().filter(|&n| n == max_letter).count();
			let normalised = self.normalise(frequency as u32, 1, column.len() as u32) as u8;

			match max_letter
			{
				b'A' => { conserved_plane.set_fg_rgb8(216, 27, 96); },
				b'T' => { conserved_plane.set_fg_rgb8(30, 136, 229); },
				b'G' => { conserved_plane.set_fg_rgb8(255, 193, 7); },
				b'C' => { conserved_plane.set_fg_rgb8(0, 77, 64); },
				b'-' => { conserved_plane.set_fg_rgb8(254, 97, 0); },
				_ => {},
			};

			let mut modified_letter = *max_letter as char;
			if normalised != 4
			{
				modified_letter.make_ascii_lowercase();
			}

			conserved_plane.putstr(modified_letter.encode_utf8(&mut tmp));

			for at in 0..normalised
			{
				match max_letter
				{
					b'A' => { histogram_plane.set_bg_rgb8(216 + (at * 10), 27 + (at * 10), 96 + (at * 10)); },
					b'T' => { histogram_plane.set_bg_rgb8(30 + (at * 10), 136 + (at * 10), 229); },
					b'G' => { histogram_plane.set_bg_rgb8(255, 193 + (at * 10), 7 + (at * 10)); },
					b'C' => { histogram_plane.set_bg_rgb8(0 + (at * 10), 77 + (at * 10), 64 + (at * 10)); },
					b'-' => { histogram_plane.set_bg_rgb8(254, 97 + (at * 10), 0 + (at * 10)); },
					_ => {},
				};

				histogram_plane.cursor_move_yx(5 - at as u32, index)?;
				histogram_plane.putstr(" ")?;
			}

			index += 1;
		}

		let mut at: usize = 0;
		for alignment in self.records.rows()
		{
			if self.top_id as usize == at
			{
				self.draw_id(&mut top_left_plane, at);
				self.draw_sequence(&mut top_sequences_plane, &alignment.to_owned(), at);
			}
			else
			{
				self.draw_id(&mut left_plane, at);
				self.draw_sequence(&mut sequences_plane, &alignment.to_owned(), at);
			}

			nc.render()?;

			at += 1;
		}

		unsafe
		{
			notcurses_mouse_enable(nc);
		}

		let mut input = NcInput::new_empty();
		loop
		{
			let key = notcurses_getc_nblock(&mut nc, &mut input);

			if key == 'Q' || key == 'q'
			{
				break;
			}

			let (dim_y, dim_x) = nc.term_dim_yx();
			match key
			{
				NCKEY_DOWN =>
				{
					if sequences_plane.y() + sequences_plane.rows() as i32 > dim_y as i32
					{
						left_plane.move_rel(-1, 0)?;
						sequences_plane.move_rel(-1, 0)?;
					}
				},
				NCKEY_UP =>
				{
					if sequences_plane.abs_y() != 3
					{
						left_plane.move_rel(1, 0)?;
						sequences_plane.move_rel(1, 0)?;
					}
				},
				NCKEY_BUTTON5 =>
				{
					if sequences_plane.x() + sequences_plane.cols() as i32 > dim_x as i32
					{
						sequences_plane.move_rel(0, -1)?;
						histogram_plane.move_rel(0, -1)?;
						top_sequences_plane.move_rel(0, -1)?;
						conserved_plane.move_rel(0, -1)?;
					}
				},
				NCKEY_BUTTON4 =>
				{
					if sequences_plane.abs_x() != self.max_name_length as u32 + 2
					{
						sequences_plane.move_rel(0, 1)?;
						histogram_plane.move_rel(0, 1)?;
						top_sequences_plane.move_rel(0, 1)?;
						conserved_plane.move_rel(0, 1)?;
					}
				},
				NCKEY_RIGHT =>
				{
					if sequences_plane.x() + sequences_plane.cols() as i32 > dim_x as i32
					{
						sequences_plane.move_rel(0, -1)?;
						histogram_plane.move_rel(0, -1)?;
						top_sequences_plane.move_rel(0, -1)?;
						conserved_plane.move_rel(0, -1)?;
					}
				},
				NCKEY_LEFT =>
				{
					if sequences_plane.abs_x() != self.max_name_length as u32 + 2
					{
						sequences_plane.move_rel(0, 1)?;
						histogram_plane.move_rel(0, 1)?;
						top_sequences_plane.move_rel(0, 1)?;
						conserved_plane.move_rel(0, 1)?;
					}
				},
				_ =>
				{
				},
			}

			rsleep![&mut nc, 0, 10];
		}

		nc.stop()?;
		Ok(())
	}

	fn normalise(&self, value: u32, min_value: u32, max_value: u32) -> u32
	{
		return 1 + (value - min_value) * (4 - 1) / (max_value - min_value);
	}
}

fn main()
{
	let yaml = load_yaml!("cli.yaml");
	let matches = App::from(yaml).get_matches();

	let source_file = match matches.value_of("SOURCE")
	{
		Some(source_file) => source_file,
		None => return,
	};

	let primary_name = match matches.value_of("primary")
	{
		Some(name) => name,
		None => "",
	};

	if let Some(view) = PolyRecord::new(source_file, primary_name)
	{
		view.draw();
	}
}
