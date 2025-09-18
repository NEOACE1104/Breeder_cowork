// import matchers from '@testing-library/jest-dom/matchers';
import { expect } from 'vitest';

// expect.extend(matchers);

// Mock i18next
vi.mock('react-i18next', () => ({
  useTranslation: () => {
    return {
      t: (str: string) => str,
      i18n: {
        changeLanguage: () => new Promise(() => {}),
      },
    };
  },
}));
