use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Debug)]
pub struct MmcifParser {
    /// Stores data categories as key-value pairs
    categories: HashMap<String, Category>,
    /// Categories built from single-value items rather than a loop
    single_value_categories: HashSet<String>,
}

#[derive(Debug)]
pub struct Category {
    /// Column names for this category
    pub headers: Vec<String>,
    /// Rows of data, where each row is a vector of values
    pub rows: Vec<Vec<String>>,
}

impl Default for MmcifParser {
    fn default() -> Self {
        Self::new()
    }
}

impl MmcifParser {
    pub fn new() -> Self {
        Self {
            categories: HashMap::new(),
            single_value_categories: HashSet::new(),
        }
    }

    /// Parse an mmCIF file from a given path
    pub fn parse_file(&mut self, path: &str) -> io::Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        self.parse_reader(reader)
    }

    /// Parse mmCIF data from any source implementing BufRead
    pub fn parse_reader<R: BufRead>(&mut self, reader: R) -> io::Result<()> {
        let mut current_category: Option<String> = None;
        let mut current_headers: Vec<String> = Vec::new();
        let mut in_loop = false;
        let mut loop_has_rows = false;
        // A tag whose value follows on a later line.
        let mut pending_item: Option<(String, String)> = None;
        // A semicolon-delimited text field being collected.
        let mut text_field: Option<(String, String, String)> = None;

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            // Inside a text field, every line is value text until a line
            // starting with ';'.
            if text_field.is_some() {
                if line.starts_with(';') {
                    let (category, field, value) = text_field.take().unwrap();
                    self.set_item(&category, &field, value);
                } else if let Some((_, _, value)) = text_field.as_mut() {
                    if !value.is_empty() {
                        value.push(' ');
                    }
                    value.push_str(trimmed);
                }
                continue;
            }

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // The value of a tag seen on an earlier line.
            if let Some((category, field)) = pending_item.take() {
                if line.starts_with(';') {
                    text_field = Some((category, field, trimmed[1..].trim().to_string()));
                    continue;
                }
                if !trimmed.starts_with('_') && !is_reserved_word(trimmed) {
                    self.set_item(&category, &field, first_value(trimmed).to_string());
                    continue;
                }
                // Otherwise the tag had no value; fall through to this line.
            }

            if trimmed.starts_with("data_") {
                // Start of a new data block
                continue;
            }

            if trimmed.starts_with("loop_") {
                in_loop = true;
                loop_has_rows = false;
                current_category = None;
                current_headers.clear();
                continue;
            }

            if trimmed.starts_with('_') {
                let parts: Vec<&str> = trimmed.splitn(2, '.').collect();
                if parts.len() != 2 {
                    continue;
                }

                let category_name = parts[0][1..].to_string(); // Remove leading underscore
                let mut rest = parts[1].splitn(2, char::is_whitespace);
                let field_name = rest.next().unwrap_or("").to_string();
                let value = rest.next().unwrap_or("").trim();

                if in_loop && !loop_has_rows {
                    // Still reading the column names of this loop.
                    if current_category.is_none() {
                        current_category = Some(category_name);
                    }
                    current_headers.push(field_name);
                } else {
                    // A tag after a loop's data rows ends the loop.
                    in_loop = false;
                    if value.is_empty() {
                        pending_item = Some((category_name, field_name));
                    } else {
                        self.set_item(&category_name, &field_name, first_value(value).to_string());
                    }
                }
                continue;
            }

            if in_loop && !current_headers.is_empty() {
                // Parse data rows
                if let Some(category_name) = &current_category {
                    loop_has_rows = true;
                    let values = parse_values(trimmed);
                    let category =
                        self.categories
                            .entry(category_name.clone())
                            .or_insert(Category {
                                headers: current_headers.clone(),
                                rows: Vec::new(),
                            });
                    category.rows.push(values);
                }
            }
        }

        // A text field left open at the end of the file.
        if let Some((category, field, value)) = text_field.take() {
            self.set_item(&category, &field, value);
        }

        Ok(())
    }

    /// Stores one single-value item, collecting a category's items into one row.
    fn set_item(&mut self, category_name: &str, field: &str, value: String) {
        if !self.single_value_categories.contains(category_name)
            && self.categories.contains_key(category_name)
        {
            return; // A loop of this category was read; leave it untouched.
        }
        self.single_value_categories
            .insert(category_name.to_string());
        let category = self
            .categories
            .entry(category_name.to_string())
            .or_insert_with(|| Category {
                headers: Vec::new(),
                rows: vec![Vec::new()],
            });
        match category.headers.iter().position(|name| name == field) {
            Some(index) => category.rows[0][index] = value,
            None => {
                category.headers.push(field.to_string());
                category.rows[0].push(value);
            }
        }
    }

    /// Get a reference to a category by name
    pub fn get_category(&self, name: &str) -> Option<&Category> {
        self.categories.get(name)
    }
}

/// Returns the first value of `text`, without its quotes.
fn first_value(text: &str) -> &str {
    let mut chars = text.chars();
    match chars.next() {
        Some(quote @ ('\'' | '"')) => {
            let rest = &text[quote.len_utf8()..];
            match rest.find(quote) {
                Some(end) => &rest[..end],
                None => rest,
            }
        }
        _ => text.split_whitespace().next().unwrap_or(""),
    }
}

/// True for the CIF reserved words that cannot be a value.
fn is_reserved_word(text: &str) -> bool {
    let lowercase = text.to_ascii_lowercase();
    ["data_", "loop_", "save_", "global_", "stop_"]
        .iter()
        .any(|word| lowercase.starts_with(word))
}

impl Category {
    /// Get a column of data by header name
    pub fn get_column(&self, header: &str) -> Option<Vec<&str>> {
        let index = self.headers.iter().position(|h| h == header)?;
        Some(
            self.rows
                .iter()
                .filter_map(|row| row.get(index).map(|s| s.as_str()))
                .collect(),
        )
    }

    /// Get a row as a map of header name to value
    pub fn get_row(&self, index: usize) -> Option<HashMap<&str, &str>> {
        self.rows.get(index).map(|row| {
            self.headers
                .iter()
                .enumerate()
                .filter_map(|(i, h)| row.get(i).map(|v| (h.as_str(), v.as_str())))
                .collect()
        })
    }
}

/// Parse a line of values, handling quoted strings and simple values
fn parse_values(line: &str) -> Vec<String> {
    let mut values = Vec::new();
    let mut current_value = String::new();
    let mut in_quotes = false;

    for c in line.chars() {
        match c {
            '"' => {
                in_quotes = !in_quotes;
                if !in_quotes {
                    values.push(current_value.clone());
                    current_value.clear();
                }
            }
            ' ' | '\t' if !in_quotes => {
                if !current_value.is_empty() {
                    values.push(current_value.clone());
                    current_value.clear();
                }
            }
            _ => current_value.push(c),
        }
    }

    if !current_value.is_empty() {
        values.push(current_value);
    }

    values
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_simple_mmcif() {
        let data = r#"data_test
_entry.id 1ABC
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
ATOM 1 C
ATOM 2 N
"#;

        let mut parser = MmcifParser::new();
        parser.parse_reader(Cursor::new(data)).unwrap();

        let atom_site = parser.get_category("atom_site").unwrap();
        assert_eq!(atom_site.headers.len(), 3);
        assert_eq!(atom_site.rows.len(), 2);
        assert_eq!(
            atom_site.get_column("group_PDB").unwrap(),
            vec!["ATOM", "ATOM"]
        );
    }

    #[test]
    fn test_parse_quoted_values() {
        let data = r#"_citation.title "Some quoted title with spaces"
loop_
_entity.type
_entity.description
"polymer" "First chain"
"non-polymer" "Second chain""#;

        let mut parser = MmcifParser::new();
        parser.parse_reader(Cursor::new(data)).unwrap();

        let entity = parser.get_category("entity").unwrap();
        assert_eq!(entity.rows.len(), 2);
        assert_eq!(
            entity.get_column("type").unwrap(),
            vec!["polymer", "non-polymer"]
        );
    }
}
