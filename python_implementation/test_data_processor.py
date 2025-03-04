import unittest
import os
import tempfile
import numpy as np
import pandas as pd
from data_processor import DataProcessor

class TestDataProcessor(unittest.TestCase):
    def setUp(self):
        """Set up test environment before each test"""
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        
        # Sample input data for testing
        self.sample_scale_data = {
            'U': [0.1, 0.2, 0.3],
            'F': [0.4, 0.5, 0.6],
            'Ca': [0.05, 0.1, 0.15]
        }
    
    def tearDown(self):
        """Clean up test environment after each test"""
        # Remove temporary directory and its contents
        for filename in os.listdir(self.test_dir):
            file_path = os.path.join(self.test_dir, filename)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(f"Error cleaning up {file_path}: {e}")

    def test_scale_to_vector(self):
        """Test SCALE output to vector conversion"""
        # Create a sample SCALE output file
        sample_file = os.path.join(self.test_dir, 'sample_scale.out')
        with open(sample_file, 'w') as f:
            f.write("""
            relative cutoff;
            
            U       0.1 0.2 0.3
            F       0.4 0.5 0.6
            Ca      0.05 0.1 0.15
            """)
        
        # Call method
        result = DataProcessor.scale_to_vector(sample_file)
        
        # Assertions
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 3)
        self.assertEqual(result['U'], [0.1, 0.2, 0.3])
        self.assertEqual(result['F'], [0.4, 0.5, 0.6])
        self.assertEqual(result['Ca'], [0.05, 0.1, 0.15])

    def test_vector_to_therm(self):
        """Test conversion of vector data to Thermochimica input"""
        # Prepare output file path
        output_base = os.path.join(self.test_dir, 'test_output')
        
        # Call method
        DataProcessor.vector_to_therm(
            self.sample_scale_data, 
            output_base, 
            '500:700:100', 
            '1:5:1'
        )
        
        # Check if files were generated
        for i in range(3):
            expected_file = f'{output_base}{i}.F90'
            self.assertTrue(os.path.exists(expected_file))
            
            # Optional: Add more detailed file content checks
            with open(expected_file, 'r') as f:
                content = f.read()
                self.assertIn(f'program {os.path.basename(output_base)}{i}', content)
                self.assertIn(f'dTemperature = {500 + i*100}', content)
                self.assertIn(f'dPressure = {i + 1}', content)

    def test_text_to_excel(self):
        """Test extracting data to Excel"""
        # Create a sample Thermochimica output file
        input_file = os.path.join(self.test_dir, 'thermo_output.txt')
        output_file = os.path.join(self.test_dir, 'output.csv')
        
        with open(input_file, 'w') as f:
            f.write("""
            1.5 mol MSFL
            80.2 Moles of pairs
            Temperature 600
            Pressure 3.5
            DEBUG
            """)
        
        # Call method
        DataProcessor.text_to_excel(input_file, output_file, 'ni T P')
        
        # Verify output
        df = pd.read_csv(output_file)
        self.assertIn('ni', df.columns)
        self.assertIn('T', df.columns)
        self.assertIn('P', df.columns)

    def test_get_ion_pair(self):
        """Test extracting cation and anion"""
        test_cases = [
            ('UF4', ('U', 'F4')),
            ('CaF2', ('Ca', 'F2')),
            ('LiF', ('Li', 'F'))
        ]
        
        for species, expected in test_cases:
            result = DataProcessor.get_ion_pair(species)
            self.assertEqual(result, expected)

    def test_parse_range_values(self):
        """Test parsing range specifications"""
        test_cases = [
            ('500:700:100', 3),  # Start:Stop:Step
            ('500:700', 3),      # Start:Stop (auto-step)
            ('500, 600, 700', 3) # Explicit values
        ]
        
        for range_str, expected_length in test_cases:
            result = DataProcessor._parse_range_values(range_str, 3)
            self.assertEqual(len(result), expected_length)

    def test_error_handling(self):
        """Test error scenarios"""
        # Test file not found
        with self.assertRaises(FileNotFoundError):
            DataProcessor.scale_to_vector('non_existent_file.txt')
        
        # Test invalid range specification
        with self.assertRaises(ValueError):
            DataProcessor._parse_range_values('invalid', 3)

if __name__ == '__main__':
    unittest.main()
