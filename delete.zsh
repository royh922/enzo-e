#!/bin/zsh

# Loop through all directories that match the pattern KH-HD-00[09][09]
for dir in KH-HD-00[0-9][0-9](/); do
  # Check if the directory exists to avoid errors
  if [ -d "$dir" ]; then
    # Extract the last two digits from the directory name
    num_part=${dir:8:2}

    # Check if the number is NOT divisible by 5
    if (( num_part % 5 != 0 )); then
      echo "Deleting $dir (number $num_part is not divisible by 5)"
      # Use rm -r to recursively delete the directory and its contents
      # Be cautious with this command.
      rm -r "$dir"
    else
      echo "Keeping $dir (number $num_part is divisible by 5)"
    fi
  fi
done

echo "Done."
