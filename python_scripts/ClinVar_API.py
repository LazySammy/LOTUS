import requests

def check_variant_in_clinvar(api_key, variant_data):
    # API endpoint URL
    url = 'https://submit.ncbi.nlm.nih.gov/api/v1/submissions/'

    # Set headers
    headers = {
        'Content-Type': 'application/json',
        'SP-API-KEY': api_key
    }

    # Prepare request payload
    payload = {
        'actions': [{
            'type': 'AddData',
            'targetDb': 'clinvar',
            'data': {
                'content': variant_data
            }
        }]
    }

    # Make POST request
    response = requests.post(url, json=payload, headers=headers)

    # Check if the variant is found in ClinVar
    if response.status_code == 200:
        print('Variant found in ClinVar!')
    else:
        print(response.status_code)
        print('Variant not found in ClinVar.')

# Example usage
api_key = '244f9c731ceb181d97280e17ebdf4b1ba909'
variant_data = '''
{
  "CHROM": "1",
  "POS": 1000000,
  "REF": "A",
  "ALT": "G"
}
'''

check_variant_in_clinvar(api_key, variant_data)