import { NavLink, Outlet } from 'react-router-dom';
import { useTranslation } from 'react-i18next';
import { Button } from '@/components/ui/button';
import { Leaf, BarChart, Dna, TestTube, Combine, FileText } from 'lucide-react';

const navItems = [
  { to: '/', label: 'nav.dashboard', icon: BarChart },
  { to: '/data-qc', label: 'nav.data_qc', icon: TestTube },
  { to: '/gp', label: 'nav.gp', icon: Dna },
  { to: '/gwas', label: 'nav.gwas', icon: Dna },
  { to: '/mating', label: 'nav.mating', icon: Combine },
  { to: '/reports', label: 'nav.reports', icon: FileText },
];

const LanguageSwitcher = () => {
  const { i18n } = useTranslation();

  const changeLanguage = (lng: 'en' | 'ko') => {
    i18n.changeLanguage(lng);
  };

  return (
    <div className="flex space-x-1">
      <Button
        variant={i18n.language === 'en' ? 'secondary' : 'ghost'}
        size="sm"
        onClick={() => changeLanguage('en')}
      >
        EN
      </Button>
      <Button
        variant={i18n.language === 'ko' ? 'secondary' : 'ghost'}
        size="sm"
        onClick={() => changeLanguage('ko')}
      >
        KO
      </Button>
    </div>
  );
};


const Layout = () => {
  const { t } = useTranslation();

  return (
    <div className="flex h-screen bg-gray-100 dark:bg-gray-900">
      <aside className="w-64 flex-shrink-0 bg-white dark:bg-gray-800 border-r border-gray-200 dark:border-gray-700 flex flex-col">
        <div className="h-16 flex items-center justify-center px-4 border-b border-gray-200 dark:border-gray-700">
          <Leaf className="h-8 w-8 text-green-600" />
          <h1 className="ml-2 text-xl font-bold text-gray-800 dark:text-white">CucumberGP</h1>
        </div>
        <nav className="flex-1 p-4 space-y-2">
          {navItems.map(({ to, label, icon: Icon }) => (
            <NavLink
              key={to}
              to={to}
              className={({ isActive }) =>
                `flex items-center px-4 py-2 text-gray-700 dark:text-gray-200 rounded-md hover:bg-gray-200 dark:hover:bg-gray-700 ${
                  isActive ? 'bg-gray-200 dark:bg-gray-700 font-semibold' : ''
                }`
              }
            >
              <Icon className="h-5 w-5 mr-3" />
              {t(label)}
            </NavLink>
          ))}
        </nav>
        <div className="p-4 border-t border-gray-200 dark:border-gray-700">
          <LanguageSwitcher />
        </div>
      </aside>
      <main className="flex-1 p-6 overflow-auto">
        <Outlet />
      </main>
    </div>
  );
};

export default Layout;
