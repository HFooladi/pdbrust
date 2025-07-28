use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Debug)]
pub struct MmcifParser {
    /// Stores data categories as key-value pairs
    categories: HashMap<String, Category>,
}

#[derive(Debug)]
pub struct Category {
    /// Column names for this category
    pub headers: Vec<String>,
    /// Rows of data, where each row is a vector of values
    pub rows: Vec<Vec<String>>,
}

impl MmcifParser {
    pub fn new() -> Self {
        Self {
            categories: HashMap::new(),
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

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            if trimmed.starts_with("data_") {
                // Start of a new data block
                continue;
            }

            if trimmed.starts_with("loop_") {
                in_loop = true;
                current_headers.clear();
                continue;
            }

            if trimmed.starts_with('_') {
                let parts: Vec<&str> = trimmed.splitn(2, '.').collect();
                if parts.len() != 2 {
                    continue;
                }

                let category_name = parts[0][1..].to_string(); // Remove leading underscore
                let field_name = parts[1].split_whitespace().next().unwrap_or("").to_string();

                if in_loop {
                    if current_category.is_none() {
                        current_category = Some(category_name.clone());
                    }
                    current_headers.push(field_name);
                } else {
                    // Handle non-loop single value items
                    let remaining = parts[1].trim();
                    let value = if remaining.starts_with('"') && remaining.ends_with('"') {
                        remaining[1..remaining.len()-1].to_string()
                    } else {
                        remaining.split_whitespace().skip(1).collect::<Vec<&str>>().join(" ")
                    };
                    let category =
                        self.categories
                            .entry(category_name.clone())
                            .or_insert(Category {
                                headers: vec![field_name.clone()],
                                rows: Vec::new(),
                            });
                    category.rows.push(vec![value]);
                }
                continue;
            }

            if in_loop && !trimmed.starts_with('_') && !current_headers.is_empty() {
                // Parse data rows
                if let Some(category_name) = &current_category {
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

        Ok(())
    }

    /// Get a reference to a category by name
    pub fn get_category(&self, name: &str) -> Option<&Category> {
        self.categories.get(name)
    }
}

impl Category {
    /// Get a column of data by header name
    pub fn get_column(&self, header: &str) -> Option<Vec<&str>> {
        let index = self.headers.iter().position(|h| h == header)?;
        Some(
            self.rows
                .iter()
                .filter_map(|row| row.get(index).map(|s| s.as_str()))
                .collect()
        )
    }

    /// Get a row as a map of header name to value
    pub fn get_row(&self, index: usize) -> Option<HashMap<&str, &str>> {
        self.rows.get(index).map(|row| {
            self.headers
                .iter()
                .enumerate()
                .filter_map(|(i, h)| {
                    row.get(i).map(|v| (h.as_str(), v.as_str()))
                })
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
