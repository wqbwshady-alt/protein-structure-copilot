from openai import OpenAI
from dotenv import load_dotenv
import os


def main():
    load_dotenv()

    api_key = os.getenv("DEEPSEEK_API_KEY")

    if not api_key:
        raise ValueError("DEEPSEEK_API_KEY not found. Please check your .env file.")

    client = OpenAI(
        api_key=api_key,
        base_url="https://api.deepseek.com"
    )

    response = client.chat.completions.create(
        model="deepseek-chat",
        messages=[
            {"role": "system", "content": "You are a structural biology assistant."},
            {"role": "user", "content": "Explain hydrophobic ligand-binding pocket in one short sentence."}
        ]
    )

    print(response.choices[0].message.content)


if __name__ == "__main__":
    main()
