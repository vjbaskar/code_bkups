#!/bin/bash

# # Example usage:
# find_substring_names "pre_name_suffix" "pre_test_suffix" "pre_demo_suffix"
find_substring_names () {
    local strings=("$@")
    local num_strings=${#strings[@]}
    local min_length=${#strings[0]}

    # Find the shortest string length
    for str in "${strings[@]}"; do
        (( ${#str} < min_length )) && min_length=${#str}
    done

    # Find common prefix
    common_prefix=""
    for ((i = 0; i < min_length; i++)); do
        char="${strings[0]:i:1}"
        for str in "${strings[@]}"; do
            [[ "${str:i:1}" != "$char" ]] && break 2
        done
        common_prefix+="$char"
    done

    # Find common suffix
    common_suffix=""
    for ((i = 0; i < min_length; i++)); do
        char="${strings[0]: -i-1:1}"
        for str in "${strings[@]}"; do
            [[ "${str: -i-1:1}" != "$char" ]] && break 2
        done
        common_suffix="$char$common_suffix"
    done

    # Extract unique substring in between
    name=()
    prefix_length=${#common_prefix}
    suffix_length=${#common_suffix}
    for str in "${strings[@]}"; do
        name+=( "${str:prefix_length:${#str}-prefix_length-suffix_length}" )
    done

    # Output results
    echo "$common_suffix $common_prefix ${name[*]}"
}
