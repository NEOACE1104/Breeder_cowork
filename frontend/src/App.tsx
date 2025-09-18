import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Layout from '@/components/layout/Layout';
import DashboardPage from '@/pages/dashboard/DashboardPage';
import DataQcPage from '@/pages/data-qc/DataQcPage';
import GpPage from '@/pages/gp/GpPage';
import GwasPage from '@/pages/gwas/GwasPage';
import MatingPage from '@/pages/mating/MatingPage';
import ReportsPage from '@/pages/reports/ReportsPage';


function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<Layout />}>
          <Route index element={<DashboardPage />} />
          <Route path="data-qc" element={<DataQcPage />} />
          <Route path="gp" element={<GpPage />} />
          <Route path="gwas" element={<GwasPage />} />
          <Route path="mating" element={<MatingPage />} />
          <Route path="reports" element={<ReportsPage />} />
        </Route>
      </Routes>
    </Router>
  );
}

export default App;
