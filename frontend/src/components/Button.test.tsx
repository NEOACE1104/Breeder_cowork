import { render, screen } from '@testing-library/react';
import { describe, it, expect } from 'vitest';
import { Button } from './ui/button';

describe('Button', () => {
  it('renders the button with children', () => {
    render(<Button>Click me</Button>);
    const buttonElement = screen.getByRole('button', { name: /click me/i });
    expect(buttonElement).toBeInTheDocument();
  });
});
